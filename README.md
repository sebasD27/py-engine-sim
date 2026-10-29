# py-engine-sim

`py-engine-sim` is a small, headless Python simulator for preliminary turbojet
and nozzle studies. It is designed for parameter sweeps, optimization, and
transparent sensitivity studies. The equations follow the same reduced-order
cycle ideas used by NASA Glenn's educational EngineSim, but this is an
independent implementation rather than a line-by-line port of the Java GUI.

## What it models

- ISA atmosphere through the lower stratosphere
- Freestream stagnation conditions
- Inlet total-pressure recovery
- Single-spool dry turbojet Brayton-cycle estimate
- Compressor work and turbine work balance
- Burner temperature and pressure loss
- Convergent and convergent-divergent nozzle flow
- Gross thrust, ram drag, pressure thrust, net thrust, fuel flow, and SFC
- JSON command-line output suitable for sweeps

The model intentionally exposes assumptions that are often unavailable for a
small commercial turbojet. It is useful for comparing candidate assumptions;
it is not a JetCat performance map and it cannot replace manufacturer data,
static testing, high-speed testing, or flight safety analysis.

## Install and run

```bash
python -m pip install -e .
py-engine-sim --preset p1000 --altitude-m 3000 --mach 1.05
```

The command prints JSON. Use `--pretty` for readable output.

```bash
py-engine-sim --preset p1000 --altitude-m 0 --mach 0 --pretty
```

## Python API

```python
from enginesim.presets import p1000_like
from enginesim.turbojet import FlightCondition, evaluate_turbojet

engine = p1000_like()
result = evaluate_turbojet(
    engine,
    FlightCondition(altitude_m=3000.0, mach=1.05),
    throttle=1.0,
)

print(result.net_thrust_n, result.fuel_flow_kg_s)
```

## Parameter calibration

The `p1000_like()` preset separates published inputs from assumptions. The
published values are used as calibration targets, not as a claim that every
internal station variable is known. A defensible workflow is:

1. Match sea-level static thrust, airflow, and fuel flow as far as the
   available data permits.
2. Sweep inlet recovery, turbine temperature, component efficiencies, and
   nozzle throat area.
3. Run altitude/Mach sweeps and report a range, not one falsely precise value.
4. Replace assumptions with manufacturer maps or test data as they become
   available.

## Units and limitations

The public API uses SI units: metres, seconds, pascals, kelvin, kilograms,
newtons, and kilograms per second. The cycle uses constant specific heats and
simple component efficiencies. It does not model a real compressor map,
surge, rotating stall, combustor dynamics, ECU control logic, inlet distortion,
3-D shocks, separated flow, or aeroelastic effects.

## Sources

See [NOTICE.md](NOTICE.md) for NASA source links and attribution. NASA's
original pages explain that EngineSim accepts altitude, Mach, throttle, inlet,
compressor, burner, turbine, and nozzle inputs, while the nozzle tool uses
throat area, exit-area ratio, chamber total conditions, and free-stream
pressure.

