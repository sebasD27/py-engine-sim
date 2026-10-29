from enginesim.atmosphere import isa


def test_sea_level_standard_atmosphere():
    state = isa(0.0)
    assert abs(state.temperature_r - 518.67) < 1e-8
    assert abs(state.pressure_psf - 2116.2166) < 0.01
    assert abs(state.density_slug_ft3 - 0.0023769) < 1e-6


def test_atmosphere_decreases_with_altitude():
    low = isa(0.0)
    high = isa(3000.0)
    assert high.pressure_psf < low.pressure_psf
    assert high.density_slug_ft3 < low.density_slug_ft3
