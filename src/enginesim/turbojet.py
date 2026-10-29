"""Reduced-order, single-spool dry turbojet model."""

from dataclasses import dataclass, asdict
import math
from typing import Any

from .atmosphere import GAMMA_AIR, R_AIR, isa, stagnation_pressure, stagnation_temperature
from .nozzle import NozzleResult, solve_nozzle


@dataclass(frozen=True)
class FlightCondition:
    altitude_m: float = 0.0
    mach: float = 0.0


@dataclass(frozen=True)
class TurbojetConfig:
    name: str = "generic-turbojet"
    compressor_pressure_ratio: float = 4.0
    design_airflow_kg_s: float = 1.8
    design_static_thrust_n: float | None = None
    design_fuel_flow_kg_s: float | None = None
    engine_mass_kg: float | None = None
    diameter_m: float | None = None
    length_m: float | None = None
    inlet_pressure_recovery: float = 0.96
    compressor_efficiency: float = 0.78
    burner_efficiency: float = 0.98
    burner_pressure_ratio: float = 0.95
    turbine_efficiency: float = 0.84
    mechanical_efficiency: float = 0.98
    turbine_inlet_temperature_k: float = 1300.0
    idle_turbine_inlet_temperature_k: float = 850.0
    nozzle_efficiency: float = 0.97
    nozzle_throat_area_m2: float = 0.015
    nozzle_exit_area_ratio: float = 1.0
    gamma_hot: float = 1.30
    gas_constant_hot_j_kgk: float = 287.0
    fuel_heating_value_j_kg: float = 43.0e6
    reference_pressure_pa: float = 101325.0
    reference_temperature_k: float = 288.15
    reference_airflow_kg_s: float = 1.8
    airflow_exponent: float = 0.85
    compressor_face_mach: float = 0.20
    max_throttle: float = 1.0
    min_throttle: float = 0.20


@dataclass(frozen=True)
class TurbojetResult:
    name: str
    altitude_m: float
    mach: float
    ambient_temperature_k: float
    ambient_pressure_pa: float
    ambient_density_kg_m3: float
    flight_speed_m_s: float
    inlet_total_pressure_pa: float
    inlet_total_temperature_k: float
    compressor_exit_pressure_pa: float
    compressor_exit_temperature_k: float
    burner_exit_pressure_pa: float
    turbine_exit_pressure_pa: float
    turbine_exit_temperature_k: float
    airflow_kg_s: float
    fuel_flow_kg_s: float
    fuel_to_air_ratio: float
    gross_thrust_n: float
    ram_drag_n: float
    pressure_thrust_n: float
    net_thrust_n: float
    specific_fuel_consumption_kg_n_h: float
    nozzle: NozzleResult

    def as_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["nozzle"] = asdict(self.nozzle)
        return data


def evaluate_turbojet(config: TurbojetConfig, flight: FlightCondition, throttle: float = 1.0) -> TurbojetResult:
    """Evaluate a reduced-order dry turbojet at one flight condition."""
    if not 0.0 <= throttle <= 1.0:
        raise ValueError("throttle must be between 0 and 1")
    if flight.mach < 0.0:
        raise ValueError("Mach number cannot be negative")
    if not 0.0 < config.inlet_pressure_recovery <= 1.0:
        raise ValueError("inlet pressure recovery must be in (0, 1]")
    atm = isa(flight.altitude_m)
    gamma = atm.gamma
    speed = flight.mach * atm.speed_of_sound_m_s
    t_t0 = stagnation_temperature(atm.temperature_k, flight.mach, gamma)
    p_t0 = stagnation_pressure(atm.pressure_pa, flight.mach, gamma)
    p_t2 = p_t0 * config.inlet_pressure_recovery
    t_t2 = t_t0

    effective_throttle = config.min_throttle + throttle * (config.max_throttle - config.min_throttle)
    pressure_scale = (p_t2 / config.reference_pressure_pa) ** config.airflow_exponent
    temperature_scale = math.sqrt(config.reference_temperature_k / t_t2)
    airflow = config.reference_airflow_kg_s * pressure_scale * temperature_scale * (0.55 + 0.45 * effective_throttle)

    cpa = 1004.5
    cpg = 1150.0
    rcp = (config.compressor_pressure_ratio - 1.0) * effective_throttle + 1.0
    t_t3s = t_t2 * rcp ** ((gamma - 1.0) / gamma)
    t_t3 = t_t2 + (t_t3s - t_t2) / config.compressor_efficiency
    p_t3 = p_t2 * rcp

    t_t4 = config.idle_turbine_inlet_temperature_k + effective_throttle * (
        config.turbine_inlet_temperature_k - config.idle_turbine_inlet_temperature_k
    )
    p_t4 = p_t3 * config.burner_pressure_ratio
    available_heat_per_kg_fuel = config.burner_efficiency * config.fuel_heating_value_j_kg - cpg * t_t4
    required_heat_per_kg_air = max(0.0, cpa * (t_t4 - t_t3))
    fuel_to_air = required_heat_per_kg_air / max(1.0, available_heat_per_kg_fuel)
    fuel_flow = airflow * fuel_to_air

    compressor_power = airflow * cpa * (t_t3 - t_t2)
    gas_flow = airflow + fuel_flow
    turbine_drop = compressor_power / max(1.0, gas_flow * cpg * config.mechanical_efficiency)
    turbine_drop = min(turbine_drop, 0.95 * t_t4)
    t_t5 = t_t4 - turbine_drop
    t_t5s = max(1.0, t_t4 - turbine_drop / config.turbine_efficiency)
    turbine_pressure_ratio = (t_t5s / t_t4) ** (config.gamma_hot / (config.gamma_hot - 1.0))
    p_t5 = max(100.0, p_t4 * turbine_pressure_ratio)

    nozzle = solve_nozzle(
        total_pressure_pa=p_t5,
        total_temperature_k=t_t5,
        ambient_pressure_pa=atm.pressure_pa,
        throat_area_m2=config.nozzle_throat_area_m2,
        exit_area_ratio=config.nozzle_exit_area_ratio,
        gamma=config.gamma_hot,
        gas_constant_j_kgk=config.gas_constant_hot_j_kgk,
        nozzle_efficiency=config.nozzle_efficiency,
        incoming_mass_flow_kg_s=airflow,
        flight_speed_m_s=speed,
    )
    gross = nozzle.gross_thrust_n * (gas_flow / max(nozzle.mass_flow_kg_s, 1e-12))
    ram_drag = airflow * speed
    net = gross - ram_drag
    sfc = fuel_flow * 3600.0 / net if net > 1e-9 else math.inf
    return TurbojetResult(
        name=config.name,
        altitude_m=flight.altitude_m,
        mach=flight.mach,
        ambient_temperature_k=atm.temperature_k,
        ambient_pressure_pa=atm.pressure_pa,
        ambient_density_kg_m3=atm.density_kg_m3,
        flight_speed_m_s=speed,
        inlet_total_pressure_pa=p_t2,
        inlet_total_temperature_k=t_t2,
        compressor_exit_pressure_pa=p_t3,
        compressor_exit_temperature_k=t_t3,
        burner_exit_pressure_pa=p_t4,
        turbine_exit_pressure_pa=p_t5,
        turbine_exit_temperature_k=t_t5,
        airflow_kg_s=airflow,
        fuel_flow_kg_s=fuel_flow,
        fuel_to_air_ratio=fuel_to_air,
        gross_thrust_n=gross,
        ram_drag_n=ram_drag,
        pressure_thrust_n=nozzle.pressure_thrust_n,
        net_thrust_n=net,
        specific_fuel_consumption_kg_n_h=sfc,
        nozzle=nozzle,
    )
