"""Small altitude/Mach sweep suitable as a starting point for optimization."""

import csv

from enginesim.presets import p1000_like
from enginesim.turbojet import FlightCondition, evaluate_turbojet


def main() -> None:
    engine = p1000_like()
    with open("p1000_sweep.csv", "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["altitude_m", "mach", "net_thrust_n", "fuel_flow_kg_s", "airflow_kg_s", "ram_drag_n"])
        for altitude in range(0, 5001, 1000):
            for step in range(0, 12):
                mach = step / 10.0
                result = evaluate_turbojet(engine, FlightCondition(altitude, mach))
                writer.writerow([
                    altitude,
                    mach,
                    result.net_thrust_n,
                    result.fuel_flow_kg_s,
                    result.airflow_kg_s,
                    result.ram_drag_n,
                ])


if __name__ == "__main__":
    main()
