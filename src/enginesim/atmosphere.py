"""International Standard Atmosphere helpers.

The model covers the troposphere and lower isothermal stratosphere, which is
adequate for the preliminary small-aircraft studies this package targets.
"""

from dataclasses import dataclass
import math


G0 = 9.80665
R_AIR = 287.05287
GAMMA_AIR = 1.4
T_SEA_LEVEL = 288.15
P_SEA_LEVEL = 101325.0
LAPSE = -0.0065
TROP_PAUSE_M = 11000.0


@dataclass(frozen=True)
class Atmosphere:
    altitude_m: float
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    speed_of_sound_m_s: float
    gamma: float = GAMMA_AIR
    gas_constant_j_kgk: float = R_AIR


def isa(altitude_m: float) -> Atmosphere:
    """Return a standard-atmosphere state for geometric altitude in metres."""
    if altitude_m < -5000.0 or altitude_m > 20000.0:
        raise ValueError("isa() supports altitudes from -5000 m through 20000 m")

    h = float(altitude_m)
    if h <= TROP_PAUSE_M:
        temperature = T_SEA_LEVEL + LAPSE * h
        pressure = P_SEA_LEVEL * (temperature / T_SEA_LEVEL) ** (
            -G0 / (LAPSE * R_AIR)
        )
    else:
        temperature = T_SEA_LEVEL + LAPSE * TROP_PAUSE_M
        pressure_tropopause = P_SEA_LEVEL * (temperature / T_SEA_LEVEL) ** (
            -G0 / (LAPSE * R_AIR)
        )
        pressure = pressure_tropopause * math.exp(
            -G0 * (h - TROP_PAUSE_M) / (R_AIR * temperature)
        )
    density = pressure / (R_AIR * temperature)
    sound = math.sqrt(GAMMA_AIR * R_AIR * temperature)
    return Atmosphere(h, temperature, pressure, density, sound)


def stagnation_temperature(static_temperature_k: float, mach: float, gamma: float = GAMMA_AIR) -> float:
    return static_temperature_k * (1.0 + 0.5 * (gamma - 1.0) * mach**2)


def stagnation_pressure(static_pressure_pa: float, mach: float, gamma: float = GAMMA_AIR) -> float:
    return static_pressure_pa * (1.0 + 0.5 * (gamma - 1.0) * mach**2) ** (gamma / (gamma - 1.0))
