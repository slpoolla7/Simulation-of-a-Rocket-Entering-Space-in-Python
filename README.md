# 🚀 Moon Mission Simulation Engine (V2.0 - Reworked)

A high-fidelity multi-phase aerospace physics engine written in Python. Moving beyond simple sub-orbital trajectories, this version tracks an entire flight profile from Earth launch to a controlled lunar landing, modelled after the historic **Apollo Saturn V** mission specifications.

This repository/branch represents **Phase 2 (Reworked)** of the simulator. It implements multi-stage mass jettison, a dual-body gravitational environment, and precise orbital insertion telemetry tracking.

---

## 🔄 What's New in Version 2.0
* **Multi-Phase Architecture:** Simulates 6 distinct, sequential mission phases from pad-clearance to lunar touchdown.
* **Dual-Body Gravity Environment:** Introduces combined Earth and Moon gravity fields using the Newton inverse-square law, mapping the transition through the Moon's Sphere of Influence (SOI).
* **Multi-Stage Mass Jettison:** Dynamically drops staging structures to calculate accurate instantaneous mass changes.
* **Telemetry Enhancements:** Tracks continuous G-load metrics and calculates precise Tsiolkovsky $\Delta v$ budgets per stage.

---

## 🌌 The 6 Flight Profiles Modelled
1. **S-IC Launch:** Stage 1 booster ignition through Earth's dense atmospheric layers.
2. **S-II Ascent:** Stage 2 burn pushing the vehicle to Low Earth Orbit (LEO) conditions ($185\text{ km}$, $\approx 7,797\text{ m/s}$).
3. **TLI Burn (S-IVB):** Translunar Injection accelerating the craft to $\approx 11,050\text{ m/s}$.
4. **Translunar Coast:** A free-flight path navigating the gravitational pull of both bodies over short time steps.
5. **LOI Burn:** Lunar Orbit Insertion deceleration into Low Lunar Orbit (LLO) at $110\text{ km}$.
6. **Lunar Descent:** Final powered brake deployment down to the Moon's surface.

---

## 🌟 Advanced Physics & Engineering Models
* **Variable Gravity Profile:** Continuously shifting Earth/Moon gravitational models using active planetary radii calculations.
* **Exponential Atmospheric Drag:** Simulates Earth-only air resistance boundaries using an $8.5\text{ km}$ atmosphere scale-height model with a hard cutoff boundary at $120\text{ km}$.
* **Orbital Mechanics Integration:** Utilizes realistic target values ($3,487\text{ m/s}$ TLI / $1,066\text{ m/s}$ LOI) calibrated to precise two-body celestial physics.

---

## 🛠️ Tech Stack & Prerequisites
* **Language:** Python 3.14+ (Leveraging advanced typing like `__future__.annotations`)
* **Core Modules Used:** `dataclasses`, `math` (Zero external installation required for core calculations)

---

## 🚀 How to Run the Simulation
1. Download or switch to this reworked repository.
2. Open your terminal and run the main entry-point script:
   ```bash
   python moon_simulation.py
   ```

---

## 📊 Telemetry Output Stride Matrix
The flight data structure tracks metrics frame-by-frame across all phases: 
`Phase | Time (s) | Altitude (m) | Velocity (m/s) | Acceleration (m/s²) | Fuel (kg) | Drag (N) | G-Force`
