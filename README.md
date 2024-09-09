# SPH

Smoothed Particle Hydrodynamics Simulation and Control System

## Introduction

Smoothed Particle Hydrodynamics (SPH) is a mesh-free particle-based method for fluid simulation.



## Installation

- `git clone [repo]`

## Requirements

- C++ 20
- OpenGL 3.3+
- AVX2 for Intrinsics

## Controls

- `Space` - Start/Pause
- `S` - One Step Forward
- `ESC` - Quit
- `R` - Reset Simulation
- `A` - Apply New Configuraion and Restart Simulation
- `T` - Restart Simulation
- `C` - Reset Viewpoint
- `Mouse Wheel` - Zoom In / Zoom Out
- `Mouse Drag` - Change Viewpoint  
- `Mouse Pan`  - Change Position

## Features

- SPH Solver
- Fluid Animation

## Architecture

- Classes
  - SPH System
  - Particle
  - Simulator
  - Rendering

### SPH System

### Particle

### Simulator

### Rendering

## References

- [Particle-Based Fluid Simulation for Interactive Applications](https://matthias-research.github.io/pages/publications/sca03.pdf)


- [A Survey on SPH Methods in Computer Graphics](https://animation.rwth-aachen.de/media/papers/77/2022-CGF-STAR_SPH.pdf)
- [sph-tutorial](https://sph-tutorial.physics-simulation.org/)
- [SPH Fluid Simulation in Python](https://www.youtube.com/watch?v=-0m05gzk8nk)
- [Real time simulation and control of Newtonian fluids using the Navier-Stokes equations](https://users.cg.tuwien.ac.at/zsolnai/gfx/fluid_control_msc_thesis/)