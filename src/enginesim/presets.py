"""Reference configurations used for demonstrations.

The P1000 preset is deliberately named ``p1000_like``: public JetCat data
does not expose the complete compressor/turbine/nozzle maps required for a
validated off-design model.
"""

from .turbojet import TurbojetConfig


def p1000_like() -> TurbojetConfig:
    return TurbojetConfig(
        name="P1000-like (screening model)",
        compressor_pressure_ratio=4.0,
        design_airflow_kg_s=1.8,
        design_static_thrust_n=1100.0,
        design_fuel_flow_kg_s=2.9 * 0.80 / 60.0,
        engine_mass_kg=11.0,
        diameter_m=0.234,
        length_m=0.444,
        inlet_pressure_recovery=0.96,
        compressor_efficiency=0.78,
        burner_efficiency=0.98,
        burner_pressure_ratio=0.95,
        turbine_efficiency=0.84,
        mechanical_efficiency=0.98,
        turbine_inlet_temperature_k=1300.0,
        idle_turbine_inlet_temperature_k=850.0,
        nozzle_efficiency=0.97,
        nozzle_throat_area_m2=0.015,
        nozzle_exit_area_ratio=1.0,
    )
