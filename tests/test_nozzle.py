from enginesim.nozzle import solve_nozzle


def test_nozzle_chokes_and_produces_positive_thrust():
    result = solve_nozzle(
        total_pressure_pa=250000.0,
        total_temperature_k=1000.0,
        ambient_pressure_pa=101325.0,
        throat_area_m2=0.01,
    )
    assert result.choked
    assert result.exit_mach == 1.0
    assert result.mass_flow_kg_s > 0.0
    assert result.gross_thrust_n > 0.0


def test_nozzle_expansion_increases_exit_mach():
    result = solve_nozzle(
        total_pressure_pa=500000.0,
        total_temperature_k=1200.0,
        ambient_pressure_pa=101325.0,
        throat_area_m2=0.01,
        exit_area_ratio=2.0,
    )
    assert result.choked
    assert result.exit_mach > 1.0
