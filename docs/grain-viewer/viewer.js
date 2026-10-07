import * as THREE from 'three';
import { TrackballControls } from './vendor/TrackballControls.js';

const DV_RAW = 0.05; // |processed / raw - 1| that raises a flag (closed raw surfaces only)
const DV_PUB = 0.2; // |processed / published - 1| that raises a flag
const REVIEW_KEY = 'grain-viewer.review.v1';
const REVIEWER_KEY = 'grain-viewer.reviewer';
const CACHE_LIMIT = 16;
const VIEW_LABELS = { raw: 'Original STL', proc: 'FEM mesh surface' };

const $ = (sel) => document.querySelector(sel);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

const state = {
  index: null,
  grains: [],
  byId: new Map(),
  list: [],
  current: null,
  view: 'overlay',
  edges: false,
  flat: false,
  axes: true,
  host: 'all',
  cohort: false,
  flagged: false,
  review: '',
  query: '',
  sort: 'id',
  reviews: {},
  token: 0,
};

// ---------------------------------------------------------------------------
// storage (per-viewer convenience only; always optional)
// ---------------------------------------------------------------------------
function storeGet(key, fallback) {
  try {
    const v = localStorage.getItem(key);
    return v == null ? fallback : JSON.parse(v);
  } catch { return fallback; }
}
function storeSet(key, value) {
  try { localStorage.setItem(key, JSON.stringify(value)); } catch { /* private mode etc. */ }
}

// ---------------------------------------------------------------------------
// formatting
// ---------------------------------------------------------------------------
const nf = (v, d = 0) => (v == null || Number.isNaN(v) ? '–' : v.toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d }));
const pct = (v, d = 1) => {
  if (v == null) return '–';
  const s = (v * 100).toFixed(d);
  return Number(s) === 0 ? `${(0).toFixed(d)}%` : `${v > 0 ? '+' : ''}${s}%`;
};
const small = (v) => (v != null && v < 1 ? `${+v.toPrecision(2)}` : nf(v, 1)); // keeps sliver edges visible
const evd = (volNm3) => Math.cbrt((6 * volNm3) / Math.PI);
const um3 = (nm3) => (nm3 == null ? '–' : (nm3 / 1e9).toPrecision(4));
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
const short = (h) => (h ? `${h.slice(0, 12)}…` : '–');
function lengthLabel(nm) {
  if (nm >= 1000) return `${+(nm / 1000).toPrecision(3)} µm`;
  return `${+nm.toPrecision(3)} nm`;
}

// ---------------------------------------------------------------------------
// flags
// ---------------------------------------------------------------------------
// Flags concern the processed mesh and call for review; notes describe the published input.
function flagsFor(g) {
  const out = [];
  const p = g.proc;
  if (!p) return [g.proc_note || 'No processed mesh.'];
  if (g.proc_note) out.push(g.proc_note);
  if (p.volume_rel_diff_vs_raw != null && Math.abs(p.volume_rel_diff_vs_raw) > DV_RAW) {
    out.push(`Processed volume differs from the raw STL by ${pct(p.volume_rel_diff_vs_raw)}.`);
  }
  if (p.volume_rel_diff_vs_published != null && Math.abs(p.volume_rel_diff_vs_published) > DV_PUB) {
    out.push(`Processed volume differs from the published volume by ${pct(p.volume_rel_diff_vs_published)}.`);
  }
  if (!p.sha256_matches_inventory) out.push('Processed mesh SHA-256 does not match its inventory entry.');
  if (p.components > 1) out.push(`Processed surface has ${p.components} disconnected pieces.`);
  return out;
}

const HOSTS = {
  PLAG: { label: 'Plagioclase host', study: 'Nikolaisen et al. 2022' },
  OPX: { label: 'Orthopyroxene host', study: 'Nikolaisen et al. 2022' },
  HEKLA: { label: 'Hekla 1991 basalt', study: 'Gergov et al. 2025' },
  VESUVIUS: { label: 'Vesuvius 1944 basalt', study: 'Gergov et al. 2025' },
};
const HOST_ORDER = Object.keys(HOSTS);

function notesFor(g) {
  const r = g.raw;
  const out = [];
  if (!r) return g.raw_note ? [g.raw_note] : [];
  if (!r.closed_manifold) {
    const d = [];
    if (r.boundary_edges) d.push(`${r.boundary_edges} open`);
    if (r.non_manifold_edges) d.push(`${r.non_manifold_edges} non-manifold`);
    if (r.misoriented_edges) d.push(`${r.misoriented_edges} inconsistently oriented`);
    out.push(`Published STL is not a closed, orientable surface (${d.join(', ')} edges), so its enclosed volume is undefined; it is drawn flat-shaded.`);
  }
  if (r.components > 1) out.push(`Published STL has ${r.components} disconnected pieces.`);
  return out;
}

// ---------------------------------------------------------------------------
// mesh loading
// ---------------------------------------------------------------------------
const cache = new Map();
async function fetchMesh(url) {
  if (cache.has(url)) {
    const hit = cache.get(url);
    cache.delete(url);
    cache.set(url, hit);
    return hit;
  }
  const promise = (async () => {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`${url}: HTTP ${res.status}`);
    let bytes = new Uint8Array(await res.arrayBuffer());
    if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      bytes = new Uint8Array(await new Response(stream).arrayBuffer());
    }
    return decode(bytes, url);
  })();
  cache.set(url, promise);
  promise.catch(() => cache.delete(url));
  while (cache.size > CACHE_LIMIT) {
    const [oldest, p] = cache.entries().next().value;
    cache.delete(oldest);
    p.then((geo) => geo.dispose(), () => {});
  }
  return promise;
}

function decode(bytes, url) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  if (String.fromCharCode(...bytes.subarray(0, 4)) !== 'GVM1') throw new Error(`${url}: not a GVM1 mesh`);
  const nV = dv.getUint32(4, true);
  const nT = dv.getUint32(8, true);
  const small = dv.getUint32(12, true) & 1;
  const lo = [0, 1, 2].map((i) => dv.getFloat32(16 + 4 * i, true));
  const step = [0, 1, 2].map((i) => dv.getFloat32(28 + 4 * i, true));
  let off = 40;
  const q = new Uint16Array(bytes.slice(off, off + nV * 6).buffer);
  off += nV * 6;
  off += (4 - (off % 4)) % 4;
  const idxBytes = bytes.slice(off, off + nT * 3 * (small ? 2 : 4)).buffer;
  const index = small ? new Uint16Array(idxBytes) : new Uint32Array(idxBytes);
  const pos = new Float32Array(nV * 3);
  for (let i = 0; i < pos.length; i++) pos[i] = lo[i % 3] + q[i] * step[i % 3];
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  geo.setIndex(new THREE.BufferAttribute(index, 1));
  geo.computeVertexNormals();
  geo.computeBoundingBox();
  geo.computeBoundingSphere();
  return geo;
}

// ---------------------------------------------------------------------------
// three.js scene
// ---------------------------------------------------------------------------
const canvas = $('#canvas');
const viewport = $('#viewport');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

const scene = new THREE.Scene();
const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.01, 1000);
camera.up.set(0, 0, 1);
scene.add(camera);
scene.add(new THREE.HemisphereLight(0xffffff, 0x8a8a8a, 1.4));
const keyLight = new THREE.DirectionalLight(0xffffff, 1.8);
keyLight.position.set(1, 1.5, 2);
camera.add(keyLight);
const fillLight = new THREE.DirectionalLight(0xffffff, 0.5);
fillLight.position.set(-2, -1, 0.5);
camera.add(fillLight);

const controls = new TrackballControls(camera, canvas);
controls.rotateSpeed = 3.2;
controls.zoomSpeed = 1.4;
controls.panSpeed = 0.9;
controls.staticMoving = true;
controls.addEventListener('change', () => requestRender());

const mats = {
  rawSolid: new THREE.MeshStandardMaterial({ roughness: 0.55, metalness: 0.05, side: THREE.DoubleSide }),
  rawGhost: new THREE.MeshStandardMaterial({ roughness: 0.6, metalness: 0, transparent: true, opacity: 0.3, depthWrite: false, side: THREE.DoubleSide }),
  procSolid: new THREE.MeshStandardMaterial({
    roughness: 0.5, metalness: 0.05, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
  }),
  rawEdges: new THREE.LineBasicMaterial({ transparent: true, opacity: 0.45 }),
  procEdges: new THREE.LineBasicMaterial({ transparent: true, opacity: 0.45 }),
};
mats.rawSolid.polygonOffset = true;
mats.rawSolid.polygonOffsetFactor = 1;
mats.rawSolid.polygonOffsetUnits = 1;

const groups = { raw: new THREE.Group(), proc: new THREE.Group() };
scene.add(groups.raw, groups.proc);
let radius = 100;

// axis triad
const triScene = new THREE.Scene();
const triCam = new THREE.OrthographicCamera(-1.6, 1.6, 1.6, -1.6, 0.1, 10);
triCam.up.copy(camera.up);
const triLabels = [];
for (const [axis, color] of [['X', 0xd9534f], ['Y', 0x4caf50], ['Z', 0x3b82f6]]) {
  const dir = new THREE.Vector3(axis === 'X' ? 1 : 0, axis === 'Y' ? 1 : 0, axis === 'Z' ? 1 : 0);
  triScene.add(new THREE.ArrowHelper(dir, new THREE.Vector3(), 1, color, 0.3, 0.16));
  const c = document.createElement('canvas');
  c.width = c.height = 64;
  const ctx = c.getContext('2d');
  ctx.font = 'bold 44px system-ui, sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillStyle = `#${color.toString(16).padStart(6, '0')}`;
  ctx.fillText(axis, 32, 34);
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), depthTest: false }));
  sprite.position.copy(dir.multiplyScalar(1.38));
  sprite.scale.setScalar(0.55);
  triScene.add(sprite);
  triLabels.push(sprite);
}

function applyTheme() {
  scene.background = new THREE.Color(css('--canvas'));
  mats.rawSolid.color.set(css('--raw'));
  mats.rawGhost.color.set(css('--raw'));
  mats.procSolid.color.set(css('--proc'));
  // Edges are a darker shade of their surface so dense meshes keep their hue.
  mats.rawEdges.color.set(css('--raw')).multiplyScalar(0.45);
  mats.procEdges.color.set(css('--proc')).multiplyScalar(0.45);
  mats.rawEdges.opacity = mats.procEdges.opacity = 0.6;
  requestRender();
}

let dirty = true;
function requestRender() { dirty = true; }

function viewports() {
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (state.current && !state.current.raw) return [{ x: 0, w, h, show: ['proc'] }]; // published as a mesh only
  if (state.view === 'split') {
    const half = Math.floor(w / 2);
    return [
      { x: 0, w: half, h, show: ['raw'] },
      { x: half, w: w - half, h, show: ['proc'] },
    ];
  }
  const show = state.view === 'overlay' ? ['raw', 'proc'] : [state.view];
  return [{ x: 0, w, h, show }];
}

function resize() {
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (canvas.width !== Math.floor(w * renderer.getPixelRatio()) || canvas.height !== Math.floor(h * renderer.getPixelRatio())) {
    renderer.setSize(w, h, false);
  }
  controls.handleResize();
  requestRender();
}
new ResizeObserver(resize).observe(viewport);

function setFrustum(aspect) {
  const half = radius * 1.15 * Math.max(1, 1 / aspect); // fit the narrower dimension
  camera.left = -half * aspect;
  camera.right = half * aspect;
  camera.top = half;
  camera.bottom = -half;
  camera.updateProjectionMatrix();
}

function render() {
  const h = canvas.clientHeight;
  renderer.setScissorTest(true);
  for (const vp of viewports()) {
    groups.raw.visible = vp.show.includes('raw');
    groups.proc.visible = vp.show.includes('proc');
    const ghost = state.view === 'overlay';
    for (const m of groups.raw.children) if (m.isMesh) m.material = ghost ? mats.rawGhost : mats.rawSolid;
    setFrustum(vp.w / Math.max(vp.h, 1));
    renderer.setViewport(vp.x, 0, vp.w, vp.h);
    renderer.setScissor(vp.x, 0, vp.w, vp.h);
    renderer.setClearColor(scene.background);
    renderer.render(scene, camera);
    if (state.axes) {
      const s = Math.min(96, Math.floor(vp.h / 4));
      triCam.position.set(0, 0, 4).applyQuaternion(camera.quaternion);
      triCam.quaternion.copy(camera.quaternion);
      renderer.setViewport(vp.x + vp.w - s - 4, 30, s, s);
      renderer.setScissor(vp.x + vp.w - s - 4, 30, s, s);
      renderer.autoClear = false;
      renderer.clearDepth();
      renderer.render(triScene, triCam);
      renderer.autoClear = true;
    }
  }
  renderer.setScissorTest(false);
  updateScaleBar(h);
}

let scaleInfo = null;
function updateScaleBar(viewH) {
  const nmPerPx = (camera.top - camera.bottom) / camera.zoom / Math.max(viewH, 1);
  const target = 110 * nmPerPx;
  const pow = 10 ** Math.floor(Math.log10(target));
  const nice = [1, 2, 5, 10].map((m) => m * pow).reduce((best, v) => (Math.abs(v - target) < Math.abs(best - target) ? v : best));
  const px = nice / nmPerPx;
  scaleInfo = { px, label: lengthLabel(nice) };
  const bar = $('#scalebar');
  bar.querySelector('.bar').style.width = `${px.toFixed(1)}px`;
  bar.querySelector('span').textContent = scaleInfo.label;
}

function loop() {
  controls.update();
  if (dirty) {
    dirty = false;
    render();
  }
  requestAnimationFrame(loop);
}

function frame(resetView) {
  const box = new THREE.Box3();
  for (const g of [groups.raw, groups.proc]) for (const m of g.children) if (m.isMesh) box.expandByObject(m);
  if (box.isEmpty()) return;
  const sphere = box.getBoundingSphere(new THREE.Sphere());
  radius = Math.max(sphere.radius, 1);
  controls.target.copy(sphere.center);
  if (resetView) {
    const dir = new THREE.Vector3(1, -1.25, 0.85).normalize();
    camera.position.copy(sphere.center).addScaledVector(dir, radius * 6);
    camera.up.set(0, 0, 1);
    camera.zoom = 1;
  } else {
    const d = camera.position.clone().sub(controls.target).normalize();
    camera.position.copy(sphere.center).addScaledVector(d, radius * 6);
  }
  camera.near = radius * 0.5;
  camera.far = radius * 12;
  camera.lookAt(sphere.center);
  camera.updateProjectionMatrix();
  controls.update();
  requestRender();
}

function edgeLines(geo, mat) {
  if (!geo.userData.wire) geo.userData.wire = new THREE.WireframeGeometry(geo);
  return new THREE.LineSegments(geo.userData.wire, mat);
}

function rebuildGroups(geos) {
  for (const key of ['raw', 'proc']) {
    groups[key].clear();
    const geo = geos[key];
    if (!geo) continue;
    groups[key].add(new THREE.Mesh(geo, key === 'raw' ? mats.rawSolid : mats.procSolid));
    if (state.edges) groups[key].add(edgeLines(geo, key === 'raw' ? mats.rawEdges : mats.procEdges));
  }
}
let currentGeos = {};

// Vertex normals are meaningless on non-orientable surfaces, so those raw meshes are always flat.
function applyShading() {
  const rawFlat = state.flat || (state.current?.raw ? !state.current.raw.closed_manifold : false);
  for (const [m, flat] of [[mats.rawSolid, rawFlat], [mats.rawGhost, rawFlat], [mats.procSolid, state.flat]]) {
    if (m.flatShading !== flat) {
      m.flatShading = flat;
      m.needsUpdate = true;
    }
  }
  requestRender();
}

// ---------------------------------------------------------------------------
// grain selection
// ---------------------------------------------------------------------------
async function showGrain(id, { resetView = true, push = true } = {}) {
  const g = state.byId.get(id);
  if (!g) return;
  state.current = g;
  const token = ++state.token;
  $('#title').textContent = g.id;
  document.title = `${g.id} · Magnetite Grain Viewer`;
  renderInfo(g);
  markSelected();
  if (push) writeHash();
  const status = $('#status');
  status.hidden = false;
  status.textContent = 'Loading meshes…';
  try {
    const base = 'data/';
    const [raw, proc] = await Promise.all([
      g.raw ? fetchMesh(base + g.raw.file) : Promise.resolve(null),
      g.proc ? fetchMesh(base + g.proc.file) : Promise.resolve(null),
    ]);
    if (token !== state.token) return;
    currentGeos = { raw, proc };
    for (const b of document.querySelectorAll('#views button')) b.disabled = !g.raw && b.dataset.view !== 'proc';
    rebuildGroups(currentGeos);
    applyShading();
    updateLabels();
    frame(resetView);
    status.hidden = true;
    prefetchNeighbours();
  } catch (err) {
    if (token !== state.token) return;
    status.textContent = `Could not load ${g.id}: ${err.message}`;
  }
}

function prefetchNeighbours() {
  const i = state.list.indexOf(state.current);
  for (const j of [i + 1, i - 1]) {
    const g = state.list[j];
    if (!g) continue;
    if (g.raw) fetchMesh(`data/${g.raw.file}`).catch(() => {});
    if (g.proc) fetchMesh(`data/${g.proc.file}`).catch(() => {});
  }
}

function step(delta) {
  if (!state.list.length) return;
  let i = state.list.indexOf(state.current);
  i = i < 0 ? 0 : (i + delta + state.list.length) % state.list.length;
  showGrain(state.list[i].id, { resetView: false });
}

function updateLabels() {
  if (state.current && !state.current.raw) {
    $('#label-left').textContent = ''; $('#label-right').textContent = '';
    $('#legend').innerHTML = `<span><i style="background:var(--proc)"></i>Published mesh (converted for merrill.jl)</span>`;
    return;
  }
  const split = state.view === 'split';
  $('#label-left').textContent = split ? VIEW_LABELS.raw : '';
  $('#label-right').textContent = split ? (state.current?.proc ? VIEW_LABELS.proc : 'No processed mesh') : '';
  const legend = $('#legend');
  const items = [];
  if (!split && state.view !== 'proc') items.push(['--raw', state.view === 'overlay' ? 'Original STL (ghost)' : 'Original STL']);
  if (!split && state.view !== 'raw') items.push(['--proc', state.current?.proc ? 'FEM mesh surface' : 'No processed mesh']);
  legend.innerHTML = items.map(([v, t]) => `<span><i style="background:var(${v})"></i>${esc(t)}</span>`).join('');
}

function setView(view) {
  state.view = view;
  for (const b of document.querySelectorAll('#views button')) b.setAttribute('aria-pressed', String(b.dataset.view === view));
  updateLabels();
  controls.handleResize();
  writeHash();
  requestRender();
}

// ---------------------------------------------------------------------------
// list + filters
// ---------------------------------------------------------------------------
const reviewStatus = (id) => state.reviews[id]?.status || 'unreviewed';

function applyFilters() {
  const q = state.query.trim().toUpperCase();
  let list = state.grains.filter((g) => (state.host === 'all' || g.host === state.host)
    && (!state.cohort || g.cohort100)
    && (!state.flagged || g._flags.length)
    && (!state.review || reviewStatus(g.id) === state.review)
    && (!q || g.id.includes(q)));
  const evsd = (g) => g.published.evsd_um ?? 0;
  const dv = (g) => (g.proc ? Math.abs(g.proc.volume_rel_diff_vs_published ?? 0) : Infinity);
  const sorters = {
    id: (a, b) => (a.host === b.host ? a.id.localeCompare(b.id) : HOST_ORDER.indexOf(a.host) - HOST_ORDER.indexOf(b.host)),
    evsd: (a, b) => evsd(a) - evsd(b),
    'evsd-desc': (a, b) => evsd(b) - evsd(a),
    dv: (a, b) => dv(b) - dv(a),
  };
  list = list.sort(sorters[state.sort]);
  state.list = list;
  renderList();
}

function renderList() {
  const ul = $('#list');
  ul.innerHTML = state.list.map((g) => {
    const rs = reviewStatus(g.id);
    return `<li role="option" data-id="${g.id}" aria-selected="${g === state.current}">
      <span class="dot ${rs}" title="${rs}"></span>
      <span class="gid">${g.id}</span>
      <span class="flag" title="${esc(g._flags.join(' '))}">${g._flags.length ? '⚑' : ''}</span>
      <span class="size muted">${g.published.evsd_um ? nf(g.published.evsd_um * 1000) : '–'}</span>
    </li>`;
  }).join('');
  $('#count').textContent = `${state.list.length} of ${state.grains.length} grains`;
  markSelected();
}

function markSelected() {
  for (const li of document.querySelectorAll('#list li')) {
    const on = li.dataset.id === state.current?.id;
    li.setAttribute('aria-selected', String(on));
    if (on) li.scrollIntoView({ block: 'nearest' });
  }
}

function renderProgress() {
  const n = state.grains.length;
  const ok = state.grains.filter((g) => reviewStatus(g.id) === 'ok').length;
  const att = state.grains.filter((g) => reviewStatus(g.id) === 'attention').length;
  $('#progress').innerHTML = `
    <div>Reviewed <b>${ok + att}</b> of ${n} · ${ok} accepted · ${att} need attention</div>
    <div class="meter"><i style="width:${(100 * ok) / n}%;background:var(--ok)"></i><i style="width:${(100 * att) / n}%;background:var(--bad)"></i></div>
    <div class="actions">
      <button class="btn" id="export">Export CSV</button>
      <label class="btn">Import CSV<input type="file" id="import" accept=".csv,text/csv" hidden></label>
      <button class="btn" id="about-btn">About</button>
    </div>`;
  $('#export').onclick = exportCsv;
  $('#import').onchange = (e) => e.target.files[0] && importCsv(e.target.files[0]);
  $('#about-btn').onclick = () => $('#about').showModal();
}

// ---------------------------------------------------------------------------
// info panel
// ---------------------------------------------------------------------------
function renderInfo(g) {
  const r = g.raw;
  const p = g.proc;
  const pub = g.published;
  const rv = state.reviews[g.id] || {};
  const flags = g._flags;
  const prov = state.index.provenance;
  if (!r) { $('#info').innerHTML = meshOnlyInfo(g, p, pub, flags, rv); bindReview(g); return; }
  const commit = prov.stl2fem_commit;
  const srcUrl = commit ?`https://github.com/duserzym/stl2fem/blob/${commit}/data/${r.source.split('/').map(encodeURIComponent).join('/')}` : null;
  const dims = (e) => (e ? e.map((v) => nf(v)).join(' × ') : '–');
  const pubVol = pub.volume_um3 != null ? pub.volume_um3 * 1e9 : null;
  const row = (label, a, b) => `<tr><th>${label}</th><td>${a}</td><td>${b}</td></tr>`;

  $('#info').innerHTML = `<div class="col">
    <section>
      <div class="head"><b>${g.id}</b><span class="muted">${HOSTS[g.host]?.label ?? g.host}</span></div>
      <div class="chips">
        ${g.cohort100 ? `<span class="chip">100-grain cohort${g.cohort_status ? ` · ${esc(g.cohort_status)}` : ''}</span>` : ''}
        ${p ? `<span class="chip proc">${esc(p.strategy)}</span>` : ''}
        ${pub.domain_state ? `<span class="chip">${esc(pub.domain_state)}</span>` : ''}
      </div>
    </section>

    <section>
      <h3>Checks</h3>
      <div class="flags">
        ${flags.length ? flags.map((f) => `<div>${esc(f)}</div>`).join('') : `<div class="none">No processed-mesh flags.</div>`}
        ${notesFor(g).map((n) => `<div class="note">${esc(n)}</div>`).join('')}
      </div>
    </section>

    <section>
      <h3>Published (Nikolaisen et al. 2022)</h3>
      <table class="kv">
        <tr><th>Volume</th><td>${pub.volume_um3 != null ? pub.volume_um3.toPrecision(4) : '–'} µm³</td></tr>
        <tr><th>EVSD</th><td>${pub.evsd_um != null ? nf(pub.evsd_um * 1000, 1) : '–'} nm</td></tr>
        <tr><th>Oblate–sphere–prolate</th><td>${pub.oblate_sphere_prolate != null ? pub.oblate_sphere_prolate.toPrecision(3) : '–'}</td></tr>
      </table>
    </section>
    </div><div class="col">
    <section>
      <h3>Geometry</h3>
      <table class="kv cmp">
        <thead><tr><th></th><th><span class="swatch" style="background:var(--raw)"></span>Original</th><th><span class="swatch" style="background:var(--proc)"></span>Processed</th></tr></thead>
        ${row('Volume (µm³)', r.volume_nm3 != null ? um3(r.volume_nm3) : '<span class="muted" title="Surface not closed/orientable">undefined</span>', p ? um3(p.volume_nm3) : '–')}
        ${row('Equiv. diameter (nm)', r.volume_nm3 != null ? nf(evd(r.volume_nm3), 1) : '–', p ? nf(evd(p.volume_nm3), 1) : '–')}
        ${row('ΔV vs published', pct(r.volume_rel_diff_vs_published), p ? `<b>${pct(p.volume_rel_diff_vs_published)}</b>` : '–')}
        ${row('ΔV vs original', '', p ? pct(p.volume_rel_diff_vs_raw) : '–')}
        ${row('Surface area (µm²)', (r.area_nm2 / 1e6).toPrecision(4), p ? (p.area_nm2 / 1e6).toPrecision(4) : '–')}
        ${row('Surface triangles', nf(r.n_triangles), p ? nf(p.n_surface_triangles) : '–')}
        ${row('Edge median (nm)', nf(r.edge_nm.median, 1), p ? `${nf(p.tet_edge_nm.median, 1)}<span class="muted"> tet</span>` : '–')}
        ${row('Edge range (nm)', `${small(r.edge_nm.min)}–${nf(r.edge_nm.max, 1)}`, p ? `${small(p.tet_edge_nm.min)}–${nf(p.tet_edge_nm.max, 1)}` : '–')}
      </table>
      <table class="kv" style="margin-top:8px">
        <tr><th>Extent, original (nm)</th><td>${dims(r.extent_nm)}</td></tr>
        ${p ? `<tr><th>Extent, processed (nm)</th><td>${dims(p.extent_nm)}</td></tr>
        <tr><th>FEM nodes / tetrahedra</th><td>${nf(p.n_nodes)} / ${nf(p.n_tets)}</td></tr>` : ''}
        <tr><th>Original surface</th><td>${r.closed_manifold ? 'closed manifold' : 'open / non-manifold'}${r.components > 1 ? `, ${r.components} pieces` : ''}</td></tr>
      </table>
    </section>

    <section class="review">
      <h3>Review</h3>
      <div class="seg" role="group" aria-label="Review status">
        <button data-status="unreviewed" aria-pressed="${!rv.status || rv.status === 'unreviewed'}">Unreviewed</button>
        <button data-status="ok" aria-pressed="${rv.status === 'ok'}">Accept</button>
        <button data-status="attention" aria-pressed="${rv.status === 'attention'}">Needs attention</button>
      </div>
      <textarea id="note" placeholder="Notes on ${g.id}…">${esc(rv.note || '')}</textarea>
      <input id="reviewer" type="search" placeholder="Reviewer name (saved in this browser)" value="${esc(storeGet(REVIEWER_KEY, ''))}">
      <div class="muted" id="saved">${rv.updated ? `Saved ${new Date(rv.updated).toLocaleString()}${rv.reviewer ? ` · ${esc(rv.reviewer)}` : ''}` : 'Notes stay in this browser until exported.'}</div>
    </section>

    <section>
      <h3>Provenance</h3>
      <div class="prov">
        <div>Original: ${srcUrl ? `<a href="${srcUrl}" target="_blank" rel="noopener">${esc(r.source)}</a>` : esc(r.source)}<br><span class="mono">sha256 ${short(r.source_sha256)}</span></div>
        ${p ? `<div>Processed: ${esc(p.source)}<br><span class="mono">sha256 ${short(p.sha256)}</span> ${p.sha256_matches_inventory ? '· matches inventory' : '· <b>inventory mismatch</b>'}</div>` : ''}
        <div>Displayed coordinates: nm, origin at raw bounding-box centre (${g.centre_nm.map((v) => nf(v, 0)).join(', ')} nm in the sample frame).</div>
      </div>
    </section>
    </div>`;
  bindReview(g);
}

function reviewHTML(g, rv) {
  return `<section class="review">
      <h3>Review</h3>
      <div class="seg" role="group" aria-label="Review status">
        <button data-status="unreviewed" aria-pressed="${!rv.status || rv.status === 'unreviewed'}">Unreviewed</button>
        <button data-status="ok" aria-pressed="${rv.status === 'ok'}">Accept</button>
        <button data-status="attention" aria-pressed="${rv.status === 'attention'}">Needs attention</button>
      </div>
      <textarea id="note" placeholder="Notes on ${g.id}…">${esc(rv.note || '')}</textarea>
      <input id="reviewer" type="search" placeholder="Reviewer name (saved in this browser)" value="${esc(storeGet(REVIEWER_KEY, ''))}">
      <div class="muted" id="saved">${rv.updated ? `Saved ${new Date(rv.updated).toLocaleString()}${rv.reviewer ? ` · ${esc(rv.reviewer)}` : ''}` : 'Notes stay in this browser until exported.'}</div>
    </section>`;
}

function bindReview(g) {
  for (const b of document.querySelectorAll('.review .seg button')) b.onclick = () => saveReview(g.id, { status: b.dataset.status });
  let t;
  $('#note').oninput = (e) => { clearTimeout(t); t = setTimeout(() => saveReview(g.id, { note: e.target.value }), 400); };
  $('#reviewer').onchange = (e) => storeSet(REVIEWER_KEY, e.target.value.trim());
}

// Grains published only as tetrahedral meshes (Gergov et al. 2025): one geometry column, published metrics.
function meshOnlyInfo(g, p, pub, flags, rv) {
  const dims = (e) => (e ? e.map((v) => nf(v)).join(' × ') : '–');
  const branch = 'https://github.com/duserzym/stl2fem/tree/Gergov2025_stl2msh';
  const gp = state.index.provenance.gergov || {};
  return `<div class="col">
    <section>
      <div class="head"><b>${g.id}</b><span class="muted">${HOSTS[g.host]?.label ?? g.host}</span></div>
      <div class="chips">
        <span class="chip proc">${esc(p.strategy)}</span>
        ${pub.ground_state ? `<span class="chip">ground state ${esc(pub.ground_state)} (published)</span>` : ''}
      </div>
    </section>
    <section>
      <h3>Checks</h3>
      <div class="flags">
        ${flags.length ? flags.map((f) => `<div>${esc(f)}</div>`).join('') : '<div class="none">No mesh flags.</div>'}
        ${notesFor(g).map((n) => `<div class="note">${esc(n)}</div>`).join('')}
      </div>
    </section>
    <section>
      <h3>Published (Gergov et al. 2025)</h3>
      <table class="kv">
        <tr><th>Volume</th><td>${pub.volume_um3 != null ? pub.volume_um3.toPrecision(4) : '–'} µm³</td></tr>
        <tr><th>Equivalent sphere diameter</th><td>${pub.evsd_um != null ? nf(pub.evsd_um * 1000, 1) : '–'} nm</td></tr>
        <tr><th>Flinn ratio</th><td>${pub.flinn_ratio != null ? pub.flinn_ratio.toPrecision(3) : '–'}</td></tr>
        <tr><th>LEM states</th><td>${esc((pub.lem_states || '–').replace(/;/g, ', '))}</td></tr>
        <tr><th>Ground state</th><td>${esc(pub.ground_state || '–')}</td></tr>
      </table>
    </section>
    </div><div class="col">
    <section>
      <h3>Geometry (merrill.jl mesh)</h3>
      <table class="kv">
        <tr><th>Volume</th><td>${um3(p.volume_nm3)} µm³</td></tr>
        <tr><th>Equivalent diameter</th><td>${nf(evd(p.volume_nm3), 1)} nm</td></tr>
        <tr><th>ΔV vs published</th><td><b>${pct(p.volume_rel_diff_vs_published)}</b></td></tr>
        <tr><th>Surface area</th><td>${(p.area_nm2 / 1e6).toPrecision(4)} µm²</td></tr>
        <tr><th>Nodes / tetrahedra</th><td>${nf(p.n_nodes)} / ${nf(p.n_tets)}</td></tr>
        <tr><th>Surface triangles</th><td>${nf(p.n_surface_triangles)}</td></tr>
        <tr><th>Edge median (range)</th><td>${nf(p.tet_edge_nm.median, 1)} nm (${small(p.tet_edge_nm.min)}–${nf(p.tet_edge_nm.max, 1)})</td></tr>
        <tr><th>Worst element shape</th><td>${p.shape_min != null ? p.shape_min.toFixed(3) : '–'} <span class="muted">(1 = regular)</span></td></tr>
        <tr><th>Extent</th><td>${dims(p.extent_nm)} nm</td></tr>
      </table>
    </section>
    ${reviewHTML(g, rv)}
    <section>
      <h3>Provenance</h3>
      <div class="prov">
        <div>Published: ${esc(p.source_pat)} in <a href="https://doi.org/10.5281/zenodo.11369780" target="_blank" rel="noopener">Zenodo 11369780</a><br><span class="mono">sha256 ${short(p.source_pat_sha256)}</span></div>
        <div>merrill.jl mesh: <a href="${branch}" target="_blank" rel="noopener">${esc(p.source)}</a><br><span class="mono">sha256 ${short(p.sha256)}</span> ${p.sha256_matches_inventory ? '· matches inventory' : '· <b>inventory mismatch</b>'}</div>
        <div>Displayed coordinates: nm, origin at the mesh bounding-box centre (${g.centre_nm.map((v) => nf(v, 0)).join(', ')} nm in the sample frame).${gp.branch_commit ? ` Built from Gergov2025_stl2msh ${gp.branch_commit.slice(0, 10)}.` : ''}</div>
      </div>
    </section>
    </div>`;
}

function saveReview(id, patch) {
  const prev = state.reviews[id] || {};
  const next = { ...prev, ...patch, updated: new Date().toISOString(), reviewer: storeGet(REVIEWER_KEY, '') || prev.reviewer || '' };
  if ((!next.status || next.status === 'unreviewed') && !next.note) delete state.reviews[id];
  else state.reviews[id] = next;
  storeSet(REVIEW_KEY, state.reviews);
  if (patch.status) {
    for (const b of document.querySelectorAll('.review .seg button')) b.setAttribute('aria-pressed', String(b.dataset.status === patch.status));
    renderList();
  }
  const saved = $('#saved');
  if (saved) saved.textContent = `Saved ${new Date(next.updated).toLocaleString()}${next.reviewer ? ` · ${next.reviewer}` : ''}`;
  renderProgress();
}

// ---------------------------------------------------------------------------
// CSV export / import
// ---------------------------------------------------------------------------
function csvCell(v) {
  const s = String(v ?? '');
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
function exportCsv() {
  const head = ['grain_id', 'host', 'review_status', 'note', 'reviewer', 'updated', 'auto_flags', 'raw_closed_manifold', 'proc_strategy', 'proc_sha256'];
  const lines = [head.join(',')];
  for (const g of state.grains) {
    const rv = state.reviews[g.id] || {};
    lines.push([g.id, g.host, rv.status || 'unreviewed', rv.note || '', rv.reviewer || '', rv.updated || '', g._flags.join(' | '), g.raw ? g.raw.closed_manifold : '', g.proc?.strategy || '', g.proc?.sha256 || ''].map(csvCell).join(','));
  }
  const name = storeGet(REVIEWER_KEY, '').replace(/[^\w-]+/g, '_');
  download(new Blob([`${lines.join('\r\n')}\r\n`], { type: 'text/csv' }), `grain_review${name ? `_${name}` : ''}_${new Date().toISOString().slice(0, 10)}.csv`);
}
function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = '';
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { cell += '"'; i++; } else if (c === '"') quoted = false; else cell += c;
    } else if (c === '"') quoted = true;
    else if (c === ',') { row.push(cell); cell = ''; } else if (c === '\n' || c === '\r') {
      if (c === '\r' && text[i + 1] === '\n') i++;
      row.push(cell); rows.push(row); row = []; cell = '';
    } else cell += c;
  }
  if (cell || row.length) { row.push(cell); rows.push(row); }
  return rows;
}
async function importCsv(file) {
  const rows = parseCsv(await file.text());
  const head = rows.shift() || [];
  const col = (name) => head.indexOf(name);
  const [ci, cs, cn, cr, cu] = ['grain_id', 'review_status', 'note', 'reviewer', 'updated'].map(col);
  if (ci < 0 || cs < 0) { alert('CSV needs grain_id and review_status columns.'); return; }
  let n = 0;
  for (const r of rows) {
    const id = r[ci];
    if (!state.byId.has(id)) continue;
    const status = r[cs];
    const note = cn >= 0 ? r[cn] : '';
    if ((!status || status === 'unreviewed') && !note) continue;
    const prev = state.reviews[id];
    const updated = (cu >= 0 && r[cu]) || new Date().toISOString();
    if (prev && prev.updated && prev.updated > updated) continue; // keep the newer local entry
    state.reviews[id] = { status, note, reviewer: cr >= 0 ? r[cr] : '', updated };
    n++;
  }
  storeSet(REVIEW_KEY, state.reviews);
  applyFilters();
  renderProgress();
  if (state.current) renderInfo(state.current);
  alert(`Imported ${n} review entr${n === 1 ? 'y' : 'ies'}.`);
}

function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// ---------------------------------------------------------------------------
// screenshot with scale bar burned in
// ---------------------------------------------------------------------------
function screenshot() {
  render();
  const dpr = renderer.getPixelRatio();
  const out = document.createElement('canvas');
  out.width = canvas.width;
  out.height = canvas.height;
  const ctx = out.getContext('2d');
  ctx.drawImage(canvas, 0, 0);
  ctx.scale(dpr, dpr);
  const ink = css('--bar');
  const font = 'system-ui, -apple-system, Segoe UI, sans-serif';
  const pill = (x, y, w, h) => { ctx.fillStyle = css('--bar-bg'); ctx.beginPath(); ctx.roundRect(x, y, w, h, 6); ctx.fill(); };
  // scale bar
  const h = canvas.clientHeight;
  const { px, label } = scaleInfo;
  ctx.font = `600 12px ${font}`;
  const boxW = Math.max(px, ctx.measureText(label).width) + 20;
  pill(14, h - 14 - 38, boxW, 38);
  ctx.strokeStyle = ink;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(25, h - 14 - 32); ctx.lineTo(25, h - 14 - 26); ctx.lineTo(25 + px - 2, h - 14 - 26); ctx.lineTo(25 + px - 2, h - 14 - 32);
  ctx.stroke();
  ctx.fillStyle = ink;
  ctx.fillText(label, 24, h - 14 - 9);
  // caption
  const g = state.current;
  const caption = `${g.id} · ${state.view === 'split' ? 'original | processed' : { raw: 'original STL', proc: 'FEM mesh surface', overlay: 'overlay (original ghost)' }[state.view]}`;
  ctx.font = `600 13px ${font}`;
  pill(14, 10, ctx.measureText(caption).width + 18, 26);
  ctx.fillStyle = ink;
  ctx.fillText(caption, 23, 28);
  out.toBlob((b) => download(b, `${g.id}_${state.view}.png`), 'image/png');
}

// ---------------------------------------------------------------------------
// URL hash
// ---------------------------------------------------------------------------
function writeHash() {
  if (!state.current) return;
  const h = `#${state.current.id}${state.view !== 'overlay' ? `/${state.view}` : ''}`;
  if (location.hash !== h) history.replaceState(null, '', h);
}
function readHash() {
  const [id, view] = decodeURIComponent(location.hash.slice(1)).toUpperCase().split('/');
  return { id, view: view?.toLowerCase() };
}

// ---------------------------------------------------------------------------
// wiring
// ---------------------------------------------------------------------------
function bindUi() {
  for (const b of document.querySelectorAll('[data-host]')) {
    b.onclick = () => {
      state.host = b.dataset.host;
      for (const o of document.querySelectorAll('[data-host]')) o.setAttribute('aria-pressed', String(o === b));
      applyFilters();
    };
  }
  $('#search').oninput = (e) => {
    state.query = e.target.value;
    applyFilters();
    const exact = state.byId.get(state.query.trim().toUpperCase());
    if (exact) showGrain(exact.id);
  };
  $('#f-cohort').onchange = (e) => { state.cohort = e.target.checked; applyFilters(); };
  $('#f-flagged').onchange = (e) => { state.flagged = e.target.checked; applyFilters(); };
  $('#f-review').onchange = (e) => { state.review = e.target.value; applyFilters(); };
  $('#sort').onchange = (e) => { state.sort = e.target.value; applyFilters(); };
  $('#list').onclick = (e) => {
    const li = e.target.closest('li');
    if (!li) return;
    showGrain(li.dataset.id, { resetView: false });
    $('#sidebar').classList.remove('open');
  };
  $('#prev').onclick = () => step(-1);
  $('#next').onclick = () => step(1);
  $('#menu').onclick = () => $('#sidebar').classList.toggle('open');
  for (const b of document.querySelectorAll('#views button')) b.onclick = () => setView(b.dataset.view);
  $('#t-edges').onchange = (e) => { state.edges = e.target.checked; rebuildGroups(currentGeos); requestRender(); };
  $('#t-flat').onchange = (e) => { state.flat = e.target.checked; applyShading(); };
  $('#t-axes').onchange = (e) => { state.axes = e.target.checked; requestRender(); };
  $('#reset').onclick = () => frame(true);
  $('#shot').onclick = screenshot;
  $('#theme').onclick = () => {
    const root = document.documentElement;
    const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    root.dataset.theme = dark ? 'light' : 'dark';
    storeSet('grain-viewer.theme', root.dataset.theme);
    applyTheme();
  };
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', applyTheme);
  window.addEventListener('hashchange', () => {
    const { id, view = 'overlay' } = readHash();
    if (view !== state.view && ['raw', 'proc', 'overlay', 'split'].includes(view)) setView(view);
    if (id && id !== state.current?.id && state.byId.has(id)) showGrain(id, { push: false });
  });
  document.addEventListener('keydown', (e) => {
    if (e.target.closest('input, textarea, select') || e.metaKey || e.ctrlKey || e.altKey) return;
    const k = e.key.toLowerCase();
    if (k === 'arrowright' || k === 'j') step(1);
    else if (k === 'arrowleft' || k === 'k') step(-1);
    else if (k in { 1: 1, 2: 1, 3: 1, 4: 1 }) setView(['raw', 'proc', 'overlay', 'split'][+k - 1]);
    else if (k === 'e') $('#t-edges').click();
    else if (k === 'f') $('#t-flat').click();
    else if (k === 'r') frame(true);
    else return;
    e.preventDefault();
  });
}

async function init() {
  const theme = storeGet('grain-viewer.theme', null);
  if (theme) document.documentElement.dataset.theme = theme;
  applyTheme();
  resize();
  requestAnimationFrame(loop);
  bindUi();
  try {
    const res = await fetch('data/index.json');
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    state.index = await res.json();
  } catch (err) {
    $('#status').hidden = false;
    $('#status').textContent = `Could not load data/index.json (${err.message}). Serve this folder over HTTP.`;
    $('#title').textContent = '—';
    return;
  }
  state.reviews = storeGet(REVIEW_KEY, {}) || {};
  state.grains = state.index.grains;
  for (const g of state.grains) {
    g._flags = flagsFor(g);
    state.byId.set(g.id, g);
  }
  const prov = state.index.provenance;
  document.querySelector('#about .dvr').textContent = `${DV_RAW * 100}%`;
  document.querySelector('#about .dvp').textContent = `${DV_PUB * 100}%`;
  $('#prov').textContent = `Data built ${state.index.generated} from stl2fem ${prov.stl2fem_commit?.slice(0, 10) ?? '?'}; inventories ${Object.values(prov.inventories).map((i) => `${i.file} (sha256 ${i.sha256.slice(0, 10)}…)`).join(', ')}${prov.inventory_repo_commit ? ` at ${prov.inventory_repo_commit.slice(0, 10)}` : ''}.`;
  applyFilters();
  renderProgress();
  const { id, view } = readHash();
  if (view && ['raw', 'proc', 'overlay', 'split'].includes(view)) setView(view);
  showGrain(state.byId.has(id) ? id : state.list[0].id, { push: false });
}

init();
