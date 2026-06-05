# 🚀 Multi-Variable Rocket Launch Simulator

A high-fidelity aerospace physics simulation engine written in Python. This application models sub-orbital rocket trajectories by calculating dynamic environmental forces in real time. Developed leveraging AI-assisted workflows, this repository represents Phase 1 of an ongoing project to construct a comprehensive aerospace telemetry simulator.

## 🌟 Advanced Physics Features
- **Dynamic Gravity Modelling:** Computes variable gravitational pull continuously using the Newton inverse-square law relative to the Earth's radius 6,371 kilometres.
- **Atmospheric Density Decay:** Simulates air resistance variations using an exponential decay model based on an 8.5 kilometers atmosphere scale height.
- **Aerodynamic Drag Calculations:** Factors in cross-sectional surface area, dynamic pressure velocities, and drag coefficients to calculate total Drag Newtons.
- **Mission Evaluation Logic:** Automatically assesses mission status milestones, tracking structural performance, launchpad clearance, and Karman Line 100 kilometers boundaries.
- **Telemetry Visualisations:** Features integrated data-plotting algorithms to review flight paths side-by-side.

## 🛠️ Tech Stack & Prerequisites
- **Language:** Python 3.14
- **Core Libraries:** `dataclasses`, `math`
- **Optional Visualisation:** `matplotlib` (For generating flight trajectory graphs)

To install the plotting dependency, run:
```bash
pip install matplotlib
```

## 🚀 How to Run the Simulation

1. Clone or download this repository.
2. Execute the script in your terminal:
```bash
python rocket_simulation.py
```
3. Follow the interactive CLI prompts to configure custom inputs (Dry mass, Fuel mass, Burn rates, Exhaust velocity, etc.) or press `Enter` to use the built-in engineering defaults.

## 📊 Sample Output Format
The simulator outputs an evaluation summary alongside a data stride matrix mapping:
`Time (s) | Altitude (m) | Velocity (m/s) | Acceleration (m/s²) | Fuel (kg) | Drag (N)`

## 🗺️ Engineering Roadmap
This engine serves as the foundational core. Future updates will focus on:
- [ ] Integrating multi-stage rocket configurations (booster separation logic).
- [ ] Adding 2D/3D visual flight path trajectory plotting.
- [ ] Incorporating wind shear profiles and variable atmospheric weather systems.
