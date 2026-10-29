from enginesim.nozzle import solve_nozzle


def test_nozzle_chokes_and_produces_positive_thrust():
    result = solve_nozzle(
        total_pressure_psf=5220.0,
        total_temperature_r=1800.0,
        ambient_pressure_psf=2116.2,
        throat_area_in2=15.5,
    )
    assert result.choked
    assert result.exit_mach == 1.0
    assert result.mass_flow_lbm_s > 0.0
    assert result.gross_thrust_lbf > 0.0


def test_nozzle_expansion_increases_exit_mach():
    result = solve_nozzle(
        total_pressure_psf=10440.0,
        total_temperature_r=2160.0,
        ambient_pressure_psf=2116.2,
        throat_area_in2=15.5,
        exit_area_ratio=2.0,
    )
    assert result.choked
    assert result.exit_mach > 1.0
