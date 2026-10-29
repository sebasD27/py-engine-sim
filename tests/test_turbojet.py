from enginesim.presets import p1000_like
from enginesim.turbojet import FlightCondition, evaluate_turbojet


def test_p1000_like_returns_finite_result():
    result = evaluate_turbojet(p1000_like(), FlightCondition(0.0, 0.0))
    assert result.net_thrust_lbf > 0.0
    assert result.airflow_lbm_s > 0.0
    assert result.fuel_flow_lbm_s > 0.0
    assert result.net_thrust_lbf == result.gross_thrust_lbf - result.ram_drag_lbf


def test_ram_drag_increases_with_mach():
    engine = p1000_like()
    static = evaluate_turbojet(engine, FlightCondition(0.0, 0.0))
    fast = evaluate_turbojet(engine, FlightCondition(0.0, 0.8))
    assert fast.ram_drag_lbf > static.ram_drag_lbf
