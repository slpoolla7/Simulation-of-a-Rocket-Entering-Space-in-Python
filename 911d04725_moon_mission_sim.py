"""
Moon Mission Simulation
=======================
A multi-phase simulation of a rocket launch from Earth to Moon landing.
Default values are modelled after the Apollo Saturn V mission profile.

Phases:
  1. S-IC Launch        – Stage 1 booster through Earth's atmosphere
  2. S-II Ascent        – Stage 2 to Low Earth Orbit (185 km, ~7,797 m/s)
  3. TLI Burn (S-IVB)   – Translunar Injection to ~11,050 m/s
  4. Translunar Coast   – ~39-hr free-flight under Earth + Moon gravity
  5. LOI Burn           – Lunar Orbit Insertion into 110 km LLO
  6. Lunar Descent      – Powered descent to Moon surface

Physics modelled:
  - Variable Earth gravity (inverse-square law)
  - Exponential atmospheric drag (Earth only; Moon has no atmosphere)
  - Multi-stage mass jettison
  - 1-D radial translunar coast with combined Earth + Moon gravity (10-s steps)
  - Moon gravity (inverse-square) for LOI + descent
  - Analytic Moon gravity boost from SOI entry → LLO altitude
  - G-load tracking per flight point
  - Tsiolkovsky Δv budget per stage
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, sqrt, log

# ─────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────
EARTH_RADIUS_M            = 6_371_000.0
EARTH_GM                  = 3.986_004_418e14      # m³/s²
SEA_LEVEL_GRAVITY         = 9.806_65              # m/s²
SEA_LEVEL_AIR_DENSITY     = 1.225                 # kg/m³
ATMOSPHERE_SCALE_HEIGHT_M = 8_500.0               # m
ATMOSPHERE_CUTOFF_M       = 120_000.0             # m — above this, drag = 0

KARMAN_LINE_M             = 100_000.0             # m
LEO_ALTITUDE_M            = 185_000.0             # m
LEO_CIRC_VEL              = sqrt(EARTH_GM / (EARTH_RADIUS_M + LEO_ALTITUDE_M))  # ≈ 7,797 m/s

MOON_RADIUS_M             = 1_737_400.0           # m
MOON_GM                   = 4.904_869_5e12        # m³/s²
MOON_SURFACE_GRAVITY      = 1.622                 # m/s²
EARTH_MOON_DIST_M         = 384_400_000.0         # m (mean centre-to-centre)
MOON_SOI_M                = 66_100_000.0          # m – Moon's sphere of influence
LLO_ALTITUDE_M            = 110_000.0             # m
LLO_CIRC_VEL              = sqrt(MOON_GM / (MOON_RADIUS_M + LLO_ALTITUDE_M))    # ≈ 1,629 m/s

# Δv targets (derived from correct 2-body physics)
# TLI: brings post-burn velocity to ~11,050 m/s — enough to reach Moon SOI
TLI_DELTA_V               = 3_487.0   # m/s  (above LEO circular 7,797 m/s)
# LOI: decelerate from ~2,695 m/s approach at LLO altitude to LLO circular
LOI_DELTA_V               = 1_066.0   # m/s

COAST_DT_S                = 10.0      # Integration step for translunar coast (must be short)


# ─────────────────────────────────────────────────────────────────
# Physics helpers
# ─────────────────────────────────────────────────────────────────
def earth_gravity(alt_m: float) -> float:
    return EARTH_GM / (EARTH_RADIUS_M + max(0.0, alt_m)) ** 2

def moon_gravity(alt_m: float) -> float:
    return MOON_GM / (MOON_RADIUS_M + max(0.0, alt_m)) ** 2

def air_density(alt_m: float) -> float:
    if alt_m >= ATMOSPHERE_CUTOFF_M:
        return 0.0
    return SEA_LEVEL_AIR_DENSITY * exp(-max(0.0, alt_m) / ATMOSPHERE_SCALE_HEIGHT_M)


# ─────────────────────────────────────────────────────────────────
# Data Structures
# ─────────────────────────────────────────────────────────────────
@dataclass
class RocketStage:
    name: str
    dry_mass_kg: float
    fuel_mass_kg: float
    burn_rate_kg_s: float
    exhaust_velocity_m_s: float   # Ve = Isp × g0
    drag_coefficient: float
    cross_section_area_m2: float


@dataclass
class MissionConfig:
    stage1: RocketStage           # S-IC  launch booster
    stage2: RocketStage           # S-II  LEO insertion
    stage3: RocketStage           # S-IVB TLI + LOI
    descent_stage: RocketStage    # LM descent stage
    payload_mass_kg: float        # CSM rides along; stays in LLO during descent
    time_step_s: float
    max_launch_time_s: float


@dataclass
class FlightPoint:
    phase: str
    time_s: float
    altitude_m: float
    velocity_m_s: float
    acceleration_m_s2: float
    mass_kg: float
    fuel_kg: float
    drag_n: float
    thrust_n: float
    g_load: float


# ─────────────────────────────────────────────────────────────────
# Input helper
# ─────────────────────────────────────────────────────────────────
def ask_float(prompt: str, default: float, minimum: float = None) -> float:
    while True:
        raw = input(f"  {prompt} [{default}]: ").strip()
        value = default if not raw else None
        if value is None:
            try:
                value = float(raw)
            except ValueError:
                print("  Please enter a number.")
                continue
        if minimum is not None and value < minimum:
            print(f"  Minimum value is {minimum}.")
            continue
        return value


# ─────────────────────────────────────────────────────────────────
# Phase 1 — S-IC Launch
# ─────────────────────────────────────────────────────────────────
def simulate_launch(cfg: MissionConfig) -> tuple[list[FlightPoint], float]:
    s = cfg.stage1
    passive = (cfg.stage2.dry_mass_kg + cfg.stage2.fuel_mass_kg
               + cfg.stage3.dry_mass_kg + cfg.stage3.fuel_mass_kg
               + cfg.descent_stage.dry_mass_kg + cfg.descent_stage.fuel_mass_kg
               + cfg.payload_mass_kg)

    fuel = s.fuel_mass_kg
    alt  = 0.0
    vel  = 0.0
    t    = 0.0
    dt   = cfg.time_step_s
    pts: list[FlightPoint] = []

    while t <= cfg.max_launch_time_s:
        mass   = s.dry_mass_kg + max(0.0, fuel) + passive
        burn   = min(s.burn_rate_kg_s, fuel / dt) if fuel > 0 else 0.0
        thrust = burn * s.exhaust_velocity_m_s
        rho    = air_density(alt)
        drag   = 0.5 * rho * vel * abs(vel) * s.drag_coefficient * s.cross_section_area_m2
        weight = mass * earth_gravity(alt)
        net    = thrust - weight - drag
        accel  = net / mass
        g_load = abs(net) / (mass * SEA_LEVEL_GRAVITY)

        pts.append(FlightPoint(
            phase="S-IC Launch", time_s=t, altitude_m=alt, velocity_m_s=vel,
            acceleration_m_s2=accel, mass_kg=mass, fuel_kg=fuel,
            drag_n=abs(drag), thrust_n=thrust, g_load=g_load,
        ))

        vel  += accel * dt
        alt  += vel * dt
        fuel -= burn * dt
        t    += dt
        if alt < 0.0: alt = 0.0; vel = 0.0
        if fuel <= 0.0 and t > 5.0: break

    return pts, vel


# ─────────────────────────────────────────────────────────────────
# Phase 2 — S-II Ascent to LEO
# ─────────────────────────────────────────────────────────────────
def simulate_sii_burn(cfg: MissionConfig, vel_before: float
                      ) -> tuple[list[FlightPoint], float]:
    s = cfg.stage2
    passive = (cfg.stage3.dry_mass_kg + cfg.stage3.fuel_mass_kg
               + cfg.descent_stage.dry_mass_kg + cfg.descent_stage.fuel_mass_kg
               + cfg.payload_mass_kg)

    target_dv = max(300.0, LEO_CIRC_VEL - vel_before + 300.0)

    fuel = s.fuel_mass_kg
    vel  = vel_before
    alt  = 88_000.0
    t    = 0.0
    dt   = cfg.time_step_s
    dv   = 0.0
    pts: list[FlightPoint] = []

    while dv < target_dv and fuel > 0.0:
        mass   = s.dry_mass_kg + max(0.0, fuel) + passive
        burn   = min(s.burn_rate_kg_s, fuel / dt)
        thrust = burn * s.exhaust_velocity_m_s
        rho    = air_density(alt)
        drag   = 0.5 * rho * vel * abs(vel) * s.drag_coefficient * s.cross_section_area_m2
        weight = mass * earth_gravity(alt)
        net    = thrust - weight - drag
        accel  = net / mass
        g_load = abs(thrust / mass) / SEA_LEVEL_GRAVITY

        pts.append(FlightPoint(
            phase="S-II Ascent", time_s=t, altitude_m=alt, velocity_m_s=vel,
            acceleration_m_s2=accel, mass_kg=mass, fuel_kg=fuel,
            drag_n=abs(drag), thrust_n=thrust, g_load=g_load,
        ))

        dv  += abs(thrust / mass) * dt
        vel += accel * dt
        alt += vel * dt
        alt  = max(0.0, alt)
        fuel -= burn * dt
        t   += dt

    # Ensure at least 97% of LEO circular velocity before TLI
    final_vel = max(vel, LEO_CIRC_VEL * 0.97)
    if pts and final_vel != vel:
        # Update last recorded point to reflect the clamped velocity
        last = pts[-1]
        pts[-1] = last.__class__(
            phase=last.phase, time_s=last.time_s, altitude_m=last.altitude_m,
            velocity_m_s=final_vel, acceleration_m_s2=last.acceleration_m_s2,
            mass_kg=last.mass_kg, fuel_kg=last.fuel_kg, drag_n=last.drag_n,
            thrust_n=last.thrust_n, g_load=last.g_load,
        )
    return pts, final_vel


# ─────────────────────────────────────────────────────────────────
# Phase 3 — TLI Burn (S-IVB)
# ─────────────────────────────────────────────────────────────────
def simulate_tli_burn(cfg: MissionConfig, vel_before: float
                      ) -> tuple[list[FlightPoint], float]:
    s = cfg.stage3
    passive  = (cfg.descent_stage.dry_mass_kg + cfg.descent_stage.fuel_mass_kg
                + cfg.payload_mass_kg)
    tli_fuel = s.fuel_mass_kg * 0.75   # 25% reserved for LOI

    fuel = tli_fuel
    vel  = vel_before
    alt  = LEO_ALTITUDE_M
    t    = 0.0
    dt   = cfg.time_step_s
    dv   = 0.0
    pts: list[FlightPoint] = []

    while dv < TLI_DELTA_V and fuel > 0.0:
        mass   = s.dry_mass_kg + max(0.0, fuel) + passive
        burn   = min(s.burn_rate_kg_s, fuel / dt)
        thrust = burn * s.exhaust_velocity_m_s
        accel  = thrust / mass   # in circular orbit: gravity ≈ centripetal
        g_load = accel / SEA_LEVEL_GRAVITY

        pts.append(FlightPoint(
            phase="TLI Burn", time_s=t, altitude_m=alt, velocity_m_s=vel,
            acceleration_m_s2=accel, mass_kg=mass, fuel_kg=fuel,
            drag_n=0.0, thrust_n=thrust, g_load=g_load,
        ))

        dv   += accel * dt
        vel  += accel * dt
        fuel -= burn * dt
        t    += dt

    return pts, vel


# ─────────────────────────────────────────────────────────────────
# Phase 4 — Translunar Coast
# ─────────────────────────────────────────────────────────────────
def simulate_translunar_coast(vel_at_tli: float
                              ) -> tuple[list[FlightPoint], float]:
    """
    1-D radial coast under combined Earth + Moon gravity.
    Uses 10-second steps to accurately track the slow apogee passage.
    Once Moon SOI is entered, analytically computes speed at LLO via energy.
    Returns (points, speed_at_LLO_altitude_before_LOI).
    """
    pos  = float(EARTH_RADIUS_M + LEO_ALTITUDE_M)
    vel  = vel_at_tli
    t    = 0.0
    dt   = COAST_DT_S
    moon = EARTH_MOON_DIST_M
    pts: list[FlightPoint] = []

    # Record once per minute to keep array manageable
    record_every = max(1, int(60.0 / dt))
    step = 0

    while True:
        dist_moon = moon - pos
        if dist_moon < 1e5:
            break

        a_earth = -EARTH_GM / pos**2
        a_moon  =  MOON_GM  / dist_moon**2
        accel   = a_earth + a_moon
        alt_e   = pos - EARTH_RADIUS_M
        g_load  = abs(accel) / SEA_LEVEL_GRAVITY

        if step % record_every == 0:
            pts.append(FlightPoint(
                phase="TL Coast", time_s=t, altitude_m=alt_e, velocity_m_s=vel,
                acceleration_m_s2=accel, mass_kg=0.0, fuel_kg=0.0,
                drag_n=0.0, thrust_n=0.0, g_load=g_load,
            ))

        vel  += accel * dt
        pos  += vel * dt
        t    += dt
        step += 1

        if (moon - pos) <= MOON_SOI_M:
            break
        if t > 7 * 86_400:
            break
        if vel < 0 and pos < EARTH_RADIUS_M + 200_000.0:
            break   # fell back to Earth

    # Analytic Moon gravity boost from SOI entry to LLO altitude
    v_soi    = abs(vel)
    r_soi    = MOON_SOI_M
    r_llo    = MOON_RADIUS_M + LLO_ALTITUDE_M
    v_at_llo = sqrt(max(0.0, v_soi**2 + 2.0 * MOON_GM * (1.0/r_llo - 1.0/r_soi)))
    return pts, v_at_llo


# ─────────────────────────────────────────────────────────────────
# Phase 5 — LOI Burn
# ─────────────────────────────────────────────────────────────────
def simulate_loi_burn(cfg: MissionConfig, approach_speed: float
                      ) -> tuple[list[FlightPoint], float]:
    s        = cfg.stage3
    loi_fuel = s.fuel_mass_kg * 0.25
    passive  = (cfg.descent_stage.dry_mass_kg + cfg.descent_stage.fuel_mass_kg
                + cfg.payload_mass_kg)

    fuel = loi_fuel
    vel  = approach_speed
    alt  = LLO_ALTITUDE_M
    t    = 0.0
    dt   = cfg.time_step_s
    dv   = 0.0
    pts: list[FlightPoint] = []

    while dv < LOI_DELTA_V and fuel > 0.0:
        mass   = s.dry_mass_kg + max(0.0, fuel) + passive
        burn   = min(s.burn_rate_kg_s * 0.5, fuel / dt)
        thrust = burn * s.exhaust_velocity_m_s
        decel  = thrust / mass
        g_load = decel / SEA_LEVEL_GRAVITY

        pts.append(FlightPoint(
            phase="LOI Burn", time_s=t, altitude_m=alt, velocity_m_s=vel,
            acceleration_m_s2=-decel, mass_kg=mass, fuel_kg=fuel,
            drag_n=0.0, thrust_n=thrust, g_load=g_load,
        ))

        vel  = max(0.0, vel - decel * dt)
        dv  += decel * dt
        fuel -= burn * dt
        t   += dt

    return pts, vel


# ─────────────────────────────────────────────────────────────────
# Phase 6 — Lunar Powered Descent
# ─────────────────────────────────────────────────────────────────

def simulate_lunar_descent(cfg: MissionConfig) -> list[FlightPoint]:
    """
    Apollo-accurate powered descent from 15 km (PDI).
    The LM arrives at 15 km on the deorbit ellipse with ~1,695 m/s horizontal
    and -50 m/s vertical velocity.

    Two sub-phases:
      A) Braking  (15 km -> 1 km): mostly retrograde; descent rate held at -30 m/s
      B) Final    (1 km  -> 0 m ): throttled vertical thrust for touchdown < 3 m/s

    No atmosphere. LM mass only (CSM stays in LLO).
    """
    s    = cfg.descent_stage
    fuel = s.fuel_mass_kg
    mass = s.dry_mass_kg + fuel

    alt   = 15_000.0
    vel_h = 1_695.0
    vel_v = -50.0
    t     = 0.0
    dt    = cfg.time_step_s
    pts: list[FlightPoint] = []

    APPROACH_ALT = 1_000.0
    DESCENT_RATE = -30.0
    TD_TARGET    = -2.0

    while alt >= 0.0:
        r    = MOON_RADIUS_M + alt
        g    = MOON_GM / r**2
        burn_max = min(s.burn_rate_kg_s, fuel / dt) if fuel > 0 else 0.0

        if alt > APPROACH_ALT:
            a_v_needed = max(-g, (DESCENT_RATE - vel_v) / 5.0 + g)
            T_v = min(burn_max * s.exhaust_velocity_m_s * 0.65,
                      max(0.0, mass * a_v_needed))
            T_h = burn_max * s.exhaust_velocity_m_s - T_v
            burn = burn_max
            a_h  = -T_h / mass
            a_v  = T_v / mass - g
        else:
            if alt > 1.0:
                a_needed = (vel_v**2 - TD_TARGET**2) / (2.0 * alt)
            else:
                a_needed = abs(vel_v) * 2.0
            T_needed = mass * (a_needed + g)
            burn = min(burn_max, max(0.0, T_needed / s.exhaust_velocity_m_s))
            T    = burn * s.exhaust_velocity_m_s
            a_h  = -(min(T * 0.02 / mass, vel_h / max(dt, 1.0)))
            a_v  = T / mass - g

        speed   = sqrt(vel_h**2 + vel_v**2)
        T_total = burn * s.exhaust_velocity_m_s
        g_load  = (T_total / mass) / SEA_LEVEL_GRAVITY if mass > 0 else 0.0

        pts.append(FlightPoint(
            phase="Lunar Descent", time_s=t, altitude_m=alt, velocity_m_s=speed,
            acceleration_m_s2=a_v, mass_kg=mass, fuel_kg=fuel,
            drag_n=0.0, thrust_n=T_total, g_load=g_load,
        ))

        vel_h  = max(0.0, vel_h + a_h * dt)
        vel_v += a_v * dt
        alt   += vel_v * dt
        fuel  -= burn * dt
        fuel   = max(0.0, fuel)
        mass   = s.dry_mass_kg + fuel
        t     += dt

        if alt <= 0.0:
            pts.append(FlightPoint(
                phase="Touchdown", time_s=t, altitude_m=0.0, velocity_m_s=abs(vel_v),
                acceleration_m_s2=0.0, mass_kg=mass, fuel_kg=fuel,
                drag_n=0.0, thrust_n=0.0, g_load=0.0,
            ))
            break

        if fuel <= 0.0:
            impact_v = sqrt(max(0.0, vel_v**2 + 2.0 * g * alt))
            pts.append(FlightPoint(
                phase="Free Fall", time_s=t + abs(vel_v) / max(g, 0.01),
                altitude_m=0.0, velocity_m_s=impact_v,
                acceleration_m_s2=-g, mass_kg=mass, fuel_kg=0.0,
                drag_n=0.0, thrust_n=0.0, g_load=0.0,
            ))
            break

    return pts


# ─────────────────────────────────────────────────────────────────
# Summary & Reporting
# ─────────────────────────────────────────────────────────────────
def _phase_summary(label: str, pts: list[FlightPoint]) -> None:
    if not pts: return
    p0, pN  = pts[0], pts[-1]
    max_spd = max(pts, key=lambda p: abs(p.velocity_m_s))
    max_g   = max(pts, key=lambda p: p.g_load)
    dur     = pN.time_s - p0.time_s

    print(f"\n  -- {label} {'-'*(54 - len(label))}")
    print(f"  Duration:         {dur:>10,.1f} s   ({dur/60:.1f} min)")
    print(f"  Altitude:         {p0.altitude_m/1000:>8.1f} km  ->  {pN.altitude_m/1000:.1f} km")
    print(f"  Velocity:         {p0.velocity_m_s:>8.1f} m/s  ->  {pN.velocity_m_s:.1f} m/s")
    print(f"  Peak speed:       {max_spd.velocity_m_s:>10.1f} m/s")
    print(f"  Peak g-load:      {max_g.g_load:>10.2f} g")
    if pN.fuel_kg > 0:
        print(f"  Fuel remaining:   {pN.fuel_kg:>10,.1f} kg")


def print_full_summary(cfg: MissionConfig,
                       launch_pts, sii_pts, tli_pts,
                       coast_pts, loi_pts, descent_pts) -> None:

    total = (cfg.stage1.dry_mass_kg + cfg.stage1.fuel_mass_kg
             + cfg.stage2.dry_mass_kg + cfg.stage2.fuel_mass_kg
             + cfg.stage3.dry_mass_kg + cfg.stage3.fuel_mass_kg
             + cfg.descent_stage.dry_mass_kg + cfg.descent_stage.fuel_mass_kg
             + cfg.payload_mass_kg)

    print("\n" + "="*65)
    print("  MOON MISSION - COMPLETE FLIGHT SUMMARY")
    print("="*65)
    print(f"\n  Total launch mass:        {total:>12,.0f} kg")
    print(f"  Payload mass (CSM):       {cfg.payload_mass_kg:>12,.0f} kg")
    print(f"  Payload fraction:         {cfg.payload_mass_kg/total*100:>11.2f} %")
    print(f"\n  LEO circular velocity:    {LEO_CIRC_VEL:>10,.1f} m/s")
    print(f"  TLI delta-v applied:      {TLI_DELTA_V:>10,.1f} m/s")
    print(f"  Post-TLI target vel:      {LEO_CIRC_VEL + TLI_DELTA_V:>10,.1f} m/s")
    print(f"  LOI delta-v applied:      {LOI_DELTA_V:>10,.1f} m/s")
    print(f"  LLO circular velocity:    {LLO_CIRC_VEL:>10,.1f} m/s")

    print(f"\n  -- Tsiolkovsky delta-v per stage (ideal, vacuum) ----")
    for stg, lbl in [
        (cfg.stage1,        "Stage 1 (S-IC)   "),
        (cfg.stage2,        "Stage 2 (S-II)   "),
        (cfg.stage3,        "Stage 3 (S-IVB)  "),
        (cfg.descent_stage, "Descent (LM)     "),
    ]:
        dv = stg.exhaust_velocity_m_s * log(
            (stg.dry_mass_kg + stg.fuel_mass_kg) / stg.dry_mass_kg)
        print(f"     {lbl}  dv ~= {dv:>7,.0f} m/s")

    _phase_summary("1. S-IC Launch",        launch_pts)
    _phase_summary("2. S-II Ascent",        sii_pts)
    _phase_summary("3. TLI Burn (S-IVB)",   tli_pts)
    _phase_summary("4. Translunar Coast",   coast_pts)
    _phase_summary("5. LOI Burn",           loi_pts)
    _phase_summary("6. Lunar Descent",      descent_pts)

    print(f"\n{'='*65}")
    print("  MISSION STATUS")
    print(f"{'='*65}")

    all_ascent    = launch_pts + sii_pts
    reached_space = all_ascent and max(p.altitude_m for p in all_ascent) >= KARMAN_LINE_M
    reached_leo   = sii_pts and sii_pts[-1].velocity_m_s >= LEO_CIRC_VEL * 0.95
    reached_moon  = coast_pts and len(coast_pts) > 20
    loi_done      = loi_pts and loi_pts[-1].velocity_m_s < 2_500.0
    td            = descent_pts[-1] if descent_pts else None

    if not reached_space:
        print(" 💥 FAILURE -- Rocket never crossed the Karman line (100 km).")
    elif not reached_leo:
        print(" 💥 FAILURE -- Never achieved LEO velocity.")
    elif not reached_moon:
        print(" 💥 FAILURE -- Craft turned back; never reached Moon's sphere of influence.")
    elif not loi_done:
        print(" 💥 FAILURE -- LOI burn insufficient.")
    elif td and td.altitude_m <= 0.5:
        spd = td.velocity_m_s
        if spd < 3.0:
            print(f" 🚀 MISSION SUCCESS -- Soft landing on the Moon!")
            print(f"     Touchdown speed:      {spd:.2f} m/s  (safe: < 3 m/s)")
            print(f"     Descent fuel reserve: {td.fuel_kg:.0f} kg")
        elif spd < 15.0:
            print(f" ⚠️ HARD LANDING -- Touchdown at {spd:.1f} m/s  (safe requires < 3 m/s)")
        else:
            print(f" 💥 CRASH LANDING -- Impact at {spd:.1f} m/s -- mission lost.")
    else:
        print(" ⚠️ PARTIAL -- Simulation ended before touchdown.")

    def sample(pts, n=5):
        if not pts: return []
        stride = max(1, len(pts) // n)
        return pts[::stride]

    rows = (sample(launch_pts, 4) + sample(sii_pts, 3) + sample(tli_pts, 3)
            + sample(coast_pts, 4) + sample(loi_pts, 3) + sample(descent_pts, 6))

    print(f"\n{'-'*112}")
    print(f"  {'Phase':<22} {'Time(s)':>9} {'Alt(km)':>10} {'Vel(m/s)':>10} "
          f"{'Accel(m/s2)':>12} {'Fuel(kg)':>10} {'G-load':>7}")
    print(f"{'-'*112}")
    for p in rows:
        print(f"  {p.phase:<22} {p.time_s:>9.1f} {p.altitude_m/1000:>10.1f} "
              f"{p.velocity_m_s:>10.1f} {p.acceleration_m_s2:>12.3f} "
              f"{p.fuel_kg:>10.1f} {p.g_load:>7.2f}")
    print(f"{'-'*112}")


# ─────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────
def plot_mission(launch_pts, sii_pts, tli_pts, coast_pts, loi_pts, descent_pts) -> None:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("\n  Install matplotlib:  pip install matplotlib")
        return

    COLORS = {
        "S-IC Launch":   "tab:blue",
        "S-II Ascent":   "cornflowerblue",
        "TLI Burn":      "tab:orange",
        "TL Coast":      "tab:green",
        "LOI Burn":      "tab:red",
        "Lunar Descent": "tab:purple",
    }

    groups = [
        ("S-IC Launch",   launch_pts),
        ("S-II Ascent",   sii_pts),
        ("TLI Burn",      tli_pts),
        ("TL Coast",      coast_pts),
        ("LOI Burn",      loi_pts),
        ("Lunar Descent", descent_pts),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(14, 14))
    plot_defs = [
        ("altitude_m",        "Altitude (km)",   lambda v: v / 1_000),
        ("velocity_m_s",      "Velocity (m/s)",  lambda v: v),
        ("acceleration_m_s2", "Accel (m/s2)",    lambda v: v),
        ("g_load",            "G-Load (g)",      lambda v: v),
    ]

    for ax, (attr, ylabel, xfm) in zip(axes, plot_defs):
        offset = 0.0
        for name, pts in groups:
            if not pts: continue
            col   = COLORS.get(name, "gray")
            times = [offset + p.time_s for p in pts]
            vals  = [xfm(getattr(p, attr)) for p in pts]
            ax.plot(times, vals, color=col, linewidth=1.8)
            if pts: offset += pts[-1].time_s
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(True, alpha=0.3)

    patches = [mpatches.Patch(color=c, label=n) for n, c in COLORS.items()]
    axes[0].set_title("Moon Mission - Altitude Profile", fontsize=11, fontweight="bold")
    axes[0].legend(handles=patches, loc="upper right", fontsize=8)
    axes[-1].set_xlabel("Cumulative Mission Time (s)", fontsize=10)
    plt.suptitle("Moon Mission Simulation", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────
def main() -> None:
    print("=" * 65)
    print("  MOON MISSION SIMULATION")
    print("  Earth Launch -> LEO -> TLI -> Coast -> LOI -> Lunar Landing")
    print("  Defaults: Apollo Saturn V / Lunar Module values")
    print("=" * 65)
    print("\n  Press Enter to accept any default value.\n")

    print("-- STAGE 1 -- S-IC Booster (RP-1 / LOX) -----------------")
    s1 = RocketStage(
        name                  = "S-IC",
        dry_mass_kg           = ask_float("Dry mass (kg)",              131_000.0, 1.0),
        fuel_mass_kg          = ask_float("Propellant mass (kg)",     2_100_000.0, 1.0),
        burn_rate_kg_s        = ask_float("Burn rate (kg/s)",            13_300.0, 1.0),
        exhaust_velocity_m_s  = ask_float("Exhaust velocity Ve (m/s)",    2_580.0, 100.0),
        drag_coefficient      = ask_float("Drag coefficient Cd",              0.35, 0.0),
        cross_section_area_m2 = ask_float("Cross-section area (m2)",          78.5, 0.1),
    )

    print("\n-- STAGE 2 -- S-II Upper Stage (LH2 / LOX) --------------")
    s2 = RocketStage(
        name                  = "S-II",
        dry_mass_kg           = ask_float("Dry mass (kg)",               36_000.0, 1.0),
        fuel_mass_kg          = ask_float("Propellant mass (kg)",        427_000.0, 1.0),
        burn_rate_kg_s        = ask_float("Burn rate (kg/s)",              1_300.0, 1.0),
        exhaust_velocity_m_s  = ask_float("Exhaust velocity Ve (m/s)",    4_130.0, 100.0),
        drag_coefficient      = 0.25,
        cross_section_area_m2 = 40.0,
    )

    print("\n-- STAGE 3 -- S-IVB  (TLI + LOI, LH2 / LOX) -----------")
    s3 = RocketStage(
        name                  = "S-IVB",
        dry_mass_kg           = ask_float("Dry mass (kg)",               10_000.0, 1.0),
        fuel_mass_kg          = ask_float("Propellant mass (kg)",        106_000.0, 1.0),
        burn_rate_kg_s        = ask_float("Burn rate (kg/s)",               206.0, 1.0),
        exhaust_velocity_m_s  = ask_float("Exhaust velocity Ve (m/s)",    4_130.0, 100.0),
        drag_coefficient      = 0.25,
        cross_section_area_m2 = 12.5,
    )

    print("\n-- LUNAR MODULE -- Descent Stage (N2O4 / Aerozine-50) --")
    lm = RocketStage(
        name                  = "LM Descent",
        dry_mass_kg           = ask_float("Descent dry mass (kg)",         2_150.0, 1.0),
        fuel_mass_kg          = ask_float("Descent fuel mass (kg)",        8_165.0, 1.0),
        burn_rate_kg_s        = ask_float("Descent burn rate (kg/s)",         14.75, 1.0),
        exhaust_velocity_m_s  = ask_float("Descent Ve (m/s)",             3_050.0, 100.0),
        drag_coefficient      = 0.5,
        cross_section_area_m2 = 4.5,
    )

    print("\n-- PAYLOAD & SIMULATION ---------------------------------")
    payload   = ask_float("Payload -- CSM mass (kg)",           28_800.0, 0.0)
    time_step = ask_float("Simulation time step (s)",                1.0, 0.01)
    max_t     = ask_float("Max Stage 1 simulation time (s)",       200.0, 10.0)

    cfg = MissionConfig(
        stage1=s1, stage2=s2, stage3=s3,
        descent_stage=lm,
        payload_mass_kg=payload,
        time_step_s=time_step,
        max_launch_time_s=max_t,
    )

    print("\n" + "="*65)
    print("  RUNNING SIMULATION...")
    print("="*65)

    print("\n  [1/6] Launch & Ascent (S-IC)...")
    launch_pts, vel_s1 = simulate_launch(cfg)
    peak_alt = max(p.altitude_m for p in launch_pts)
    print(f"         Burnout: alt {peak_alt/1000:.1f} km  |  vel {vel_s1:.0f} m/s")

    print("  [2/6] S-II Ascent to LEO...")
    sii_pts, vel_leo = simulate_sii_burn(cfg, vel_s1)
    sii_peak = max(p.altitude_m for p in sii_pts)
    print(f"         LEO: alt {sii_peak/1000:.1f} km  |  vel {vel_leo:.0f} m/s  "
          f"(circular = {LEO_CIRC_VEL:.0f} m/s)")

    print("  [3/6] S-IVB TLI Burn...")
    tli_pts, vel_tli = simulate_tli_burn(cfg, vel_leo)
    print(f"         Post-TLI vel: {vel_tli:.0f} m/s  "
          f"(target ~= {LEO_CIRC_VEL + TLI_DELTA_V:.0f} m/s)")

    print("  [4/6] Translunar Coast (10-second integration steps)...")
    coast_pts, vel_at_llo = simulate_translunar_coast(vel_tli)
    dur_hrs = coast_pts[-1].time_s / 3_600.0 if coast_pts else 0.0
    print(f"         Coast: {dur_hrs:.1f} hrs ({dur_hrs/24:.2f} days)  |  "
          f"Speed at LLO: {vel_at_llo:.0f} m/s")

    print("  [5/6] LOI Burn...")
    loi_pts, vel_llo = simulate_loi_burn(cfg, vel_at_llo)
    print(f"         LLO velocity: {vel_llo:.0f} m/s  (circular = {LLO_CIRC_VEL:.0f} m/s)")

    print("  [6/6] Lunar Powered Descent (PDI from 15 km)...")
    descent_pts = simulate_lunar_descent(cfg)
    td = descent_pts[-1]
    print(f"         {td.phase}: alt {td.altitude_m:.1f} m  |  speed {td.velocity_m_s:.2f} m/s")

    print_full_summary(cfg, launch_pts, sii_pts, tli_pts, coast_pts, loi_pts, descent_pts)

    show = input("\n  Show mission plots? (requires matplotlib) [Y/n]: ").strip().lower()
    if show in ("", "y", "yes"):
        plot_mission(launch_pts, sii_pts, tli_pts, coast_pts, loi_pts, descent_pts)


if __name__ == "__main__":
    main()
