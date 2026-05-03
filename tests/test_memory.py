from stl2fem.memory import estimate_merrill_memory


def test_memory_estimate_scales_with_mesh_size():
    small = estimate_merrill_memory(100, 500)
    large = estimate_merrill_memory(200, 1000)

    assert small["estimated_memory_bytes"] > 0
    assert large["estimated_memory_bytes"] > small["estimated_memory_bytes"]
    assert large["estimated_memory_human"].endswith(("KiB", "MiB", "GiB"))

