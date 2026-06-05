from dataclasses import dataclass
from math import exp

# Physics Constants 
EARTH_RADIUS_M = 6_371_000.0
SEA_LEVEL_GRAVITY = 9.80665
SEA_LEVEL_AIR_DENSITY = 1.225
ATMOSPHERE_SCALE_HEIGHT_M = 8_500.0
KARMAN_LINE_M = 100_000.0  # Definition of where space begins


# Data Structures
@dataclass
class RocketConfig:
    dry_mass_kg: float
    fuel_mass_kg: float
    burn_rate_kg_s: float
    exhaust_velocity_m_s: float
    drag_coefficient: float
    cross_section_area_m2: float
    time_step_s: float
    max_time_s: float


@dataclass
class FlightPoint:
    time_s: float
    altitude_m: float
    velocity_m_s: float
    acceleration_m_s2: float
    mass_kg: float
    fuel_kg: float
    drag_n: float
    thrust_n: float


# Helper Functions
def ask_float(prompt: str, default: float, minimum: float = None) -> float:
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        if not raw:
            value = default
        else:
            try:
                value = float(raw)
            except ValueError:
                print("Please enter a number.")
                continue

        if minimum is not None and value < minimum:
            print(f"Please enter a value of at least {minimum}.")
            continue
        return value


def gravity_at_altitude(altitude_m: float) -> float:
    distance_from_center = EARTH_RADIUS_M + max(0.0, altitude_m)
    return SEA_LEVEL_GRAVITY * (EARTH_RADIUS_M / distance_from_center) ** 2


def air_density_at_altitude(altitude_m: float) -> float:
    return SEA_LEVEL_AIR_DENSITY * exp(-max(0.0, altitude_m) / ATMOSPHERE_SCALE_HEIGHT_M)


# Simulation Core
def simulate(config: RocketConfig) -> list[FlightPoint]:
    fuel_kg = config.fuel_mass_kg
    dry_mass_kg = config.dry_mass_kg
    altitude_m = 0.0
    velocity_m_s = 0.0
    time_s = 0.0
    points: list[FlightPoint] = []

    while time_s <= config.max_time_s:
        mass_kg = dry_mass_kg + fuel_kg
        actual_burn_rate = min(config.burn_rate_kg_s, fuel_kg / config.time_step_s)
        thrust_n = actual_burn_rate * config.exhaust_velocity_m_s

        density = air_density_at_altitude(altitude_m)
        drag_n = (
            0.5
            * density
            * velocity_m_s
            * abs(velocity_m_s)
            * config.drag_coefficient
            * config.cross_section_area_m2
        )
        weight_n = mass_kg * gravity_at_altitude(altitude_m)
        net_force_n = thrust_n - weight_n - drag_n
        acceleration_m_s2 = net_force_n / mass_kg

        points.append(
            FlightPoint(
                time_s=time_s,
                altitude_m=altitude_m,
                velocity_m_s=velocity_m_s,
                acceleration_m_s2=acceleration_m_s2,
                mass_kg=mass_kg,
                fuel_kg=fuel_kg,
                drag_n=abs(drag_n),
                thrust_n=thrust_n,
            )
        )

        if altitude_m <= 0.0 and time_s > 2.0 and fuel_kg <= 0.0:
            break

        velocity_m_s += acceleration_m_s2 * config.time_step_s
        altitude_m += velocity_m_s * config.time_step_s
        fuel_kg -= actual_burn_rate * config.time_step_s
        time_s += config.time_step_s

        if altitude_m < 0.0:
            altitude_m = 0.0
            velocity_m_s = 0.0

    return points


def print_summary(points: list[FlightPoint], config: RocketConfig) -> None:
    max_altitude_point = max(points, key=lambda point: point.altitude_m)
    max_speed_point = max(points, key=lambda point: abs(point.velocity_m_s))
    burnout = next((point for point in points if point.fuel_kg <= 0.0), None)
    
    # Mission evaluation at the end of the simulation
    final_point = points[-1]
    has_landed_or_crashed = final_point.altitude_m <= 0.1 and final_point.time_s > 5.0

    print("\nFlight summary")
    print("-" * 48)
    print(f"Lift-off mass:        {config.dry_mass_kg + config.fuel_mass_kg:,.1f} kg")
    print(f"Initial thrust:       {config.burn_rate_kg_s * config.exhaust_velocity_m_s:,.1f} N")
    print(f"Peak altitude:        {max_altitude_point.altitude_m:,.1f} m at {max_altitude_point.time_s:.1f} s")
    print(f"Peak speed:           {max_speed_point.velocity_m_s:,.1f} m/s at {max_speed_point.time_s:.1f} s")
    if burnout:
        print(f"Fuel burnout:         {burnout.time_s:.1f} s at {burnout.altitude_m:,.1f} m")
    else:
        print("Fuel burnout:         not reached during simulation")
    print(f"Final altitude:       {final_point.altitude_m:,.1f} m")
    print(f"Final velocity:       {final_point.velocity_m_s:,.1f} m/s")

    #
    Mission Evaluation Logic
    print("\nMission Status")
    print("-" * 48)
    
    # 1. Did it ever clear the launchpad?
    if max_altitude_point.altitude_m < 10.0:
        print("❌ CRITICAL FAILURE: The rocket failed to leave the launchpad (insufficient Thrust-to-Weight ratio).")
    
    # 2. Did it reach space?
    elif max_altitude_point.altitude_m >= KARMAN_LINE_M:
        if has_landed_or_crashed:
            print("🚀 SUCCESS: The rocket successfully reached space (crossed the Karman Line)!")
            print("💥 Note: The booster has impacted the ground after a successful sub-orbital flight.")
        else:
            print("🚀 SUCCESS: The rocket is currently in space or still in high-altitude flight!")
            
    # 3. Did it fly but fail to reach space?
    else:
        print("⚠️ SUB-ORBITAL TEST: The rocket launched successfully but did not reach space.")
        if has_landed_or_crashed:
            print("💥 Status: The vehicle has crashed back to Earth.")
        else:
            print("Status: The vehicle is still in mid-air at the end of the simulation timeframe.")

    print("\nSample trajectory")
    print("-" * 82)
    print(" time(s) | altitude(m) | velocity(m/s) | accel(m/s^2) | fuel(kg) | drag(N)")
    print("-" * 82)
    stride = max(1, len(points) // 12)
    for point in points[::stride]:
        print(
            f"{point.time_s:8.1f} | "
            f"{point.altitude_m:11.1f} | "
            f"{point.velocity_m_s:13.1f} | "
            f"{point.acceleration_m_s2:12.2f} | "
            f"{point.fuel_kg:8.1f} | "
            f"{point.drag_n:7.1f}"
        )


def plot_flight(points: list[FlightPoint]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nInstall matplotlib to see graphs: python3 -m pip install matplotlib")
        return

    times = [point.time_s for point in points]
    altitudes = [point.altitude_m for point in points]
    velocities = [point.velocity_m_s for point in points]
    accelerations = [point.acceleration_m_s2 for point in points]

    _, axes = plt.subplots(3, 1, sharex=True, figsize=(9, 8))
    axes[0].plot(times, altitudes, color="tab:blue")
    axes[0].set_ylabel("Altitude (m)")
    axes[0].grid(True)

    axes[1].plot(times, velocities, color="tab:green")
    axes[1].set_ylabel("Velocity (m/s)")
    axes[1].grid(True)

    axes[2].plot(times, accelerations, color="tab:red")
    axes[2].set_ylabel("Acceleration (m/s^2)")
    axes[2].set_xlabel("Time (s)")
    axes[2].grid(True)

    plt.suptitle("Rocket launch simulation")
    plt.tight_layout()
    plt.show()


def main() -> None:
    print("Rocket launch simulator")
    print("Press Enter to use the value in brackets.\n")

    config = RocketConfig(
        dry_mass_kg=ask_float("Rocket dry mass, without fuel, in kg", 12_000.0, 1.0),
        fuel_mass_kg=ask_float("Fuel mass in kg", 38_000.0, 0.0),
        burn_rate_kg_s=ask_float("Fuel burn rate in kg/s", 250.0, 0.0),
        exhaust_velocity_m_s=ask_float("Exhaust velocity in m/s", 3_000.0, 0.0),
        drag_coefficient=ask_float("Drag coefficient", 0.5, 0.0),
        cross_section_area_m2=ask_float("Cross-section area in square meters", 10.0, 0.0),
        time_step_s=ask_float("Simulation time step in seconds", 0.1, 0.001),
        max_time_s=ask_float("Maximum simulation time in seconds", 300.0, 1.0),
    )

    points = simulate(config)
    print_summary(points, config)

    show_plot = input("\nShow graphs if matplotlib is available? [Y/n]: ").strip().lower()
    if show_plot in ("", "y", "yes"):
        plot_flight(points)


if __name__ == "__main__":
    main()
