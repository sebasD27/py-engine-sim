"""One-dimensional compressible nozzle calculations."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class NozzleResult:
    choked: bool
    exit_mach: float
    exit_pressure_pa: float
    exit_temperature_k: float
    exit_velocity_m_s: float
    mass_flow_kg_s: float
    gross_thrust_n: float
    pressure_thrust_n: float
    exit_area_m2: float


def _area_ratio(mach: float, gamma: float) -> float:
    term = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (1.0 / mach) * ((2.0 / (gamma + 1.0)) * term) ** exponent


def _mach_for_area_ratio(area_ratio: float, gamma: float, supersonic: bool) -> float:
    if area_ratio < 1.0:
        raise ValueError("exit/throat area ratio must be at least one")
    if abs(area_ratio - 1.0) < 1e-12:
        return 1.0
    lo, hi = (1.0, 20.0) if supersonic else (1e-8, 1.0)
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        value = _area_ratio(mid, gamma)
        if supersonic:
            if value > area_ratio:
                lo = mid
            else:
                hi = mid
        else:
            if value > area_ratio:
                lo = mid
            else:
                hi = mid
    return 0.5 * (lo + hi)


def solve_nozzle(
    total_pressure_pa: float,
    total_temperature_k: float,
    ambient_pressure_pa: float,
    throat_area_m2: float,
    *,
    exit_area_ratio: float = 1.0,
    gamma: float = 1.30,
    gas_constant_j_kgk: float = 287.0,
    nozzle_efficiency: float = 0.97,
    incoming_mass_flow_kg_s: float | None = None,
    flight_speed_m_s: float = 0.0,
) -> NozzleResult:
    """Solve an idealized convergent or convergent-divergent nozzle.

    ``total_pressure_pa`` and ``total_temperature_k`` are nozzle-entry total
    conditions. ``incoming_mass_flow_kg_s`` is used only to report gross
    thrust minus ram drag outside this function; the returned gross thrust is
    nozzle gross thrust including pressure thrust.
    """
    if total_pressure_pa <= 0 or total_temperature_k <= 0 or ambient_pressure_pa <= 0:
        raise ValueError("pressures and temperature must be positive")
    if throat_area_m2 <= 0 or exit_area_ratio < 1.0:
        raise ValueError("throat area must be positive and exit area ratio >= 1")
    if not 0 < nozzle_efficiency <= 1:
        raise ValueError("nozzle efficiency must be in (0, 1]")

    critical_ratio = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    pressure_ratio = ambient_pressure_pa / total_pressure_pa
    choked = pressure_ratio <= critical_ratio
    mdot_factor = math.sqrt(gamma / gas_constant_j_kgk) / math.sqrt(total_temperature_k)

    if not choked:
        exit_pressure = ambient_pressure_pa
        exit_mach = _mach_from_pressure_ratio(pressure_ratio, gamma)
        exit_area = throat_area_m2
        mass_flow = (
            total_pressure_pa
            * mdot_factor
            * exit_mach
            * (1.0 + 0.5 * (gamma - 1.0) * exit_mach**2) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
            * exit_area
        )
    else:
        mass_flow = (
            total_pressure_pa
            * mdot_factor
            * math.sqrt(2.0 / (gamma + 1.0))
            * (2.0 / (gamma + 1.0)) ** (1.0 / (gamma - 1.0))
            * throat_area_m2
        )
        exit_area = throat_area_m2 * exit_area_ratio
        if exit_area_ratio > 1.0:
            exit_mach = _mach_for_area_ratio(exit_area_ratio, gamma, supersonic=True)
            pressure_ratio_exit = (1.0 + 0.5 * (gamma - 1.0) * exit_mach**2) ** (-gamma / (gamma - 1.0))
            exit_pressure = total_pressure_pa * pressure_ratio_exit
        else:
            exit_mach = 1.0
            exit_pressure = total_pressure_pa * critical_ratio

    exit_temperature = total_temperature_k / (1.0 + 0.5 * (gamma - 1.0) * exit_mach**2)
    exit_velocity = exit_mach * math.sqrt(gamma * gas_constant_j_kgk * exit_temperature)
    exit_velocity *= math.sqrt(nozzle_efficiency)
    pressure_thrust = (exit_pressure - ambient_pressure_pa) * exit_area
    gross_thrust = mass_flow * exit_velocity + pressure_thrust
    return NozzleResult(
        choked=choked,
        exit_mach=exit_mach,
        exit_pressure_pa=exit_pressure,
        exit_temperature_k=exit_temperature,
        exit_velocity_m_s=exit_velocity,
        mass_flow_kg_s=mass_flow,
        gross_thrust_n=gross_thrust,
        pressure_thrust_n=pressure_thrust,
        exit_area_m2=exit_area,
    )


def _mach_from_pressure_ratio(pressure_ratio: float, gamma: float) -> float:
    """Subsonic Mach from static-to-total pressure ratio."""
    if not 0 < pressure_ratio <= 1.0:
        raise ValueError("pressure ratio must be in (0, 1]")
    lo, hi = 0.0, 1.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        calculated = (1.0 + 0.5 * (gamma - 1.0) * mid**2) ** (-gamma / (gamma - 1.0))
        if calculated > pressure_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
