/*
    Smoothed Particle Hydrodynamics (SPH)

    Learn how to design and implement a fluid dynamics simulator from scratch in C++ using OpenGL.
        Computational Fluid Dynamics 
            Particle-based Incompressible Inviscid Laminar Fluid Flow Simulator

        Lagrangian Navier Stokes Equations

    ---------------
    Data Structures
    ---------------

        struct particle {
            v3f position;
            v3f velocity;
            v3f acceleration;
            f32 mass;
            f32 density;
            f32 pressure;
            f32 viscosity;
            f32 radius;
        };

        std::vector<Particle> particles;

    ---------------
    Algorithms
    ---------------

        SPH Algorithm
            0. Initialize Particles
            1. Compute Particle Density
            2. Compute Particle Pressure
            3. Compute Pressure Force
            4. Compute Viscous Force
            5. Compute Surface Tension Force
            6. Compute Gravity Force
            7. Update Acceleration
            8. Update Velocity
            9. Update Position

        Smoothing Kernel
            Density:   Poly6 Kernel
            Pressure:  Spikey Kernel Grad 
            Viscosity: Spikey Kernel Laplacian
            Surface Tension: 

        Neighborhood Search
            Spatial Hash Grid

        Marching Cubes
        Point Splatting
        Bubbles
        Foam

    ----------
    References
    ----------
        SPH
            1. Particle-Based Fluid Simulation for Interactive Applications: https://matthias-research.github.io/pages/publications/sca03.pdf
            2. Versatile Surface Tension and Adhesion for SPH Fluids: https://cg.informatik.uni-freiburg.de/publications/2013_SIGGRAPHASIA_surfaceTensionAdhesion.pdf
            3. GPU Fluid Simulation: https://wickedengine.net/2018/05/scalabe-gpu-fluid-simulation
            
            Smoothed Particle Hydrodynamics Fluid Simulation: https://rlguy.com/sphfluidsim/
            
            sph-tutorial: https://sph-tutorial.physics-simulation.org

            SPH Fluid Simulation in Python: https://www.youtube.com/watch?v=-0m05gzk8nk
            Smoothed-particle Hydrodynamics: https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics

            A Survey on SPH Methods in Computer Graphics: https://animation.rwth-aachen.de/media/papers/77/2022-CGF-STAR_SPH.pdf
            Particle-based Viscoelastic Fluid Simulation: https://www.ljll.fr/~frey/papers/levelsets/Clavet%20S.,%20Particle-based%20viscoelastic%20fluid%20simulation.pdf
            Unified Spray, Foam and Bubbles for Particle-Based Fluids: https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf

        Compiler Intrinsics
           x86 intrinsics list: https://learn.microsoft.com/en-us/cpp/intrinsics/x86-intrinsics-list?view=msvc-170
           Performance Optimisation of Smoothed particle Hydrodynamics Algorithms fro Multi/Many-Core Architectures: https://arxiv.org/pdf/1612.06090

        Paralell Computing
            OpenMP: https://www.openmp.org/
            Guide into OpenMP: Easy multithreading programming for C++: https://bisqwit.iki.fi/story/howto/openmp/#PrefaceImportanceOfMultithreading
            An Introduction to Parallel Computing in C++: https://www.cs.cmu.edu/afs/cs/academic/class/15210-f15/www/pasl.html#ch:race-conditions
*/
#pragma once

#include <vector>
#include <iostream>
#include <limits>

#include <immintrin.h>
//#include <omp.h>

#include "Common.h"
#include "Random.h"
#include "Particle.h"
#include "Config.h"

struct SPHSettings
{
    // Simulation
    s32 n_particles;
    f32 rest_density;
    f32 gas_constant;
    f32 viscosity;
    f32 surface_tension_constant;
    vf3 gravity;
    f32 mass;
    f32 dt;

    // Smoothing Kernel
    f32 smoothing_radius;
    f32 smoothing_radius2;
    f32 poly6;
    f32 spiky_grad;
    f32 spiky_laplacian;

    // Boundary
    f32 boundary_epsilon;
    f32 boundary_damping;
    vf3 boundary_size;
    vf3 boundary_min;
    vf3 boundary_max;
};

class SPH
{
public:
    SPH();
    SPH(const SPHSettings& s);

public:
    void Init(const SPHSettings& s);
    void Simulate(f32 dt);

public:
    void ComputeDensityPressure();
    void ComputeDensityPressureIntrinsics();

    void ComputePressureViscousSurfaceTensionForce();
    void ComputePressureViscousSurfaceTensionForceIntrinsics();

    void Integrate(f32 dt);
    void ComputeBoundaryCondition();

public:
    std::vector<Particle>& GetParticles();

private:
    SPHSettings settings;
    std::vector<Particle> fluid_particles;
    //std::vector<particle> boundary_particles; // TODO: will be introduced

private: // Spatial hash map
    void CreateNeighborTable();
    std::vector<std::vector<u32>> particle_table;
    vi3 cell(const Particle& p, f32 h);
    u32 hash(const vi3& cell);
};