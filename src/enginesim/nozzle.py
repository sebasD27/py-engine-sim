"""One-dimensional compressible nozzle calculations.

The public interface uses psf, degrees Rankine, inches squared, ft/s, lbm/s,
and lbf. The calculation itself is performed in SI units.
"""

from dataclasses import dataclass
import math

from .units import (
    FT_LBF_PER_SLUG_R_TO_J_PER_KG_K,
    FT_TO_M,
    IN_TO_M,
    K_TO_R,
    LBM_TO_KG,
    LBF_TO_N,
    PSF_TO_PA,
)


@dataclass(frozen=True)
class NozzleResult:
    choked: bool
    exit_mach: float
    exit_pressure_psf: float
    exit_temperature_r: float
    exit_velocity_ft_s: float
    mass_flow_lbm_s: float
    gross_thrust_lbf: float
    pressure_thrust_lbf: float
    exit_area_in2: float


@dataclass(frozen=True)
class _NozzleResultSI:
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
        if value > area_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve_nozzle(
    total_pressure_psf: float,
    total_temperature_r: float,
    ambient_pressure_psf: float,
    throat_area_in2: float,
    *,
    exit_area_ratio: float = 1.0,
    gamma: float = 1.30,
    gas_constant_ft_lbf_slug_r: float = 1716.59,
    nozzle_efficiency: float = 0.97,
    incoming_mass_flow_lbm_s: float | None = None,
    flight_speed_ft_s: float = 0.0,
) -> NozzleResult:
    """Solve an idealized convergent or convergent-divergent nozzle.

    Inputs and outputs use imperial aerospace units. ``incoming_mass_flow_lbm_s``
    and ``flight_speed_ft_s`` are accepted for API symmetry; nozzle gross
    thrust is calculated from nozzle flow.
    """
    if total_pressure_psf <= 0 or total_temperature_r <= 0 or ambient_pressure_psf <= 0:
        raise ValueError("pressures and temperature must be positive")
    if throat_area_in2 <= 0 or exit_area_ratio < 1.0:
        raise ValueError("throat area must be positive and exit area ratio >= 1")
    if not 0 < nozzle_efficiency <= 1:
        raise ValueError("nozzle efficiency must be in (0, 1]")

    result = _solve_nozzle_si(
        total_pressure_pa=total_pressure_psf * PSF_TO_PA,
        total_temperature_k=total_temperature_r / K_TO_R,
        ambient_pressure_pa=ambient_pressure_psf * PSF_TO_PA,
        throat_area_m2=throat_area_in2 * IN_TO_M**2,
        exit_area_ratio=exit_area_ratio,
        gamma=gamma,
        gas_constant_j_kgk=gas_constant_ft_lbf_slug_r * FT_LBF_PER_SLUG_R_TO_J_PER_KG_K,
        nozzle_efficiency=nozzle_efficiency,
    )
    return NozzleResult(
        choked=result.choked,
        exit_mach=result.exit_mach,
        exit_pressure_psf=result.exit_pressure_pa / PSF_TO_PA,
        exit_temperature_r=result.exit_temperature_k * K_TO_R,
        exit_velocity_ft_s=result.exit_velocity_m_s / FT_TO_M,
        mass_flow_lbm_s=result.mass_flow_kg_s / LBM_TO_KG,
        gross_thrust_lbf=result.gross_thrust_n / LBF_TO_N,
        pressure_thrust_lbf=result.pressure_thrust_n / LBF_TO_N,
        exit_area_in2=result.exit_area_m2 / IN_TO_M**2,
    )


def _solve_nozzle_si(
    total_pressure_pa: float,
    total_temperature_k: float,
    ambient_pressure_pa: float,
    throat_area_m2: float,
    *,
    exit_area_ratio: float,
    gamma: float,
    gas_constant_j_kgk: float,
    nozzle_efficiency: float,
) -> _NozzleResultSI:
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
            * (1.0 + 0.5 * (gamma - 1.0) * exit_mach**2) ** (
                -(gamma + 1.0) / (2.0 * (gamma - 1.0))
            )
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
            pressure_ratio_exit = (1.0 + 0.5 * (gamma - 1.0) * exit_mach**2) ** (
                -gamma / (gamma - 1.0)
            )
            exit_pressure = total_pressure_pa * pressure_ratio_exit
        else:
            exit_mach = 1.0
            exit_pressure = total_pressure_pa * critical_ratio

    exit_temperature = total_temperature_k / (1.0 + 0.5 * (gamma - 1.0) * exit_mach**2)
    exit_velocity = exit_mach * math.sqrt(gamma * gas_constant_j_kgk * exit_temperature)
    exit_velocity *= math.sqrt(nozzle_efficiency)
    pressure_thrust = (exit_pressure - ambient_pressure_pa) * exit_area
    gross_thrust = mass_flow * exit_velocity + pressure_thrust
    return _NozzleResultSI(
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
