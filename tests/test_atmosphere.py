from enginesim.atmosphere import isa


def test_sea_level_standard_atmosphere():
    state = isa(0.0)
    assert abs(state.temperature_k - 288.15) < 1e-8
    assert abs(state.pressure_pa - 101325.0) < 1e-5
    assert abs(state.density_kg_m3 - 1.225) < 0.002


def test_atmosphere_decreases_with_altitude():
    low = isa(0.0)
    high = isa(3000.0)
    assert high.pressure_pa < low.pressure_pa
    assert high.density_kg_m3 < low.density_kg_m3
