/*
    Smoothed Particle Hydrodynamics (SPH)

    Learn how to design and implement a fluid dynamics simulator from scratch in C++ using OpenGL.
        Computational Fluid Dynamics Simulator
            Particle-based Incompressible Inviscous Laminar Fluid

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
            f32 surface_tension;

            std::vector<particle*> neighbours;
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

        Smoothing Kernel Function
            Density  : Poly6 Kernel
            Pressure : Spikey Grad Kernel
            Viscosity: Spikey Laplacian Kernel
            Surface Tension : 

        Density Function
        Pressure Function
        Viscosity Function
        Surface Tension Function
        Gravity Function
        Acceleration Function
        Velocity Function
        Position Function

        Neighborhood Search
            Hash map
            Spatial Grid

        Marching Cubes
        Point Splatting
        Foam

        Time Integration (Symplectic  Euler Integration)

        

    ----------
    References
    ----------
        SPH
            1. Particle-Based Fluid Simulation for Interactive Applications: https://matthias-research.github.io/pages/publications/sca03.pdf
            2. Versatile Surface Tension and Adhesion for SPH Fluids: https://cg.informatik.uni-freiburg.de/publications/2013_SIGGRAPHASIA_surfaceTensionAdhesion.pdf
            3. GPU Fluid Simulation: https://wickedengine.net/2018/05/scalabe-gpu-fluid-simulation
            4. Smoothed Particle Hydrodynamics Fluid Simulation: https://rlguy.com/sphfluidsim/
            sph-tutorial: https://sph-tutorial.physics-simulation.org

            SPH Fluid Simulation in Python: https://www.youtube.com/watch?v=-0m05gzk8nk
            Wiki, Smoothed-particle Hydrodynamics: https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics

            A Survey on SPH Methods in Computer Graphics: https://animation.rwth-aachen.de/media/papers/77/2022-CGF-STAR_SPH.pdf
            Particle-based Viscoelastic Fluid Simulation: https://www.ljll.fr/~frey/papers/levelsets/Clavet%20S.,%20Particle-based%20viscoelastic%20fluid%20simulation.pdf
            Unified Spray, Foam and Bubbles for Particle-Based Fluids: https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf

        Paralell Computing
            OpenMP: https://www.openmp.org/
            Guide into OpenMP: Easy multithreading programming for C++: https://bisqwit.iki.fi/story/howto/openmp/#PrefaceImportanceOfMultithreading
            An Introduction to Parallel Computing in C++: https://www.cs.cmu.edu/afs/cs/academic/class/15210-f15/www/pasl.html#ch:race-conditions

*/
#pragma once

#include <vector>
#include <iostream>

#define _OPENMP_LLVM_RUNTIME
#include <omp.h>

#include "Common.h"
#include "Random.h"

const u32 SCREEN_WIDTH  = 1280;
const u32 SCREEN_HEIGHT = 720;

static const f32 PI = 3.14159265358979323846f;

struct settings
{
    // Simulation
    s32 n_particles  = 25000;
    f32 rest_density = 997.0f;
    f32 gas_constant = 1.0f;
    f32 viscosity = 0.002f;
    f32 surface_tension_constant = 0.0000001f;
    vf3 gravity = { 0.0f, -9.80665f, 0.0f };
    f32 mass = 0.01f;
    // Smoothing Kernel
    f32 smoothing_radius = 0.04f;
    f32 smoothing_radius2 = std::pow(0.04f, 2.0f);
    f32 poly6 = 315.0f / (64.0f * PI * std::pow(0.04f, 9.0f));
    f32 spiky_grad = -45.0f / (PI * std::pow(0.04f, 6.0f));
    f32 spiky_laplacian = 45.0f / (PI * std::pow(0.04f, 6.0f));
    // Boundary
    f32 boundary_epsilon = 0.0000001f;
    f32 boundary_damping = -0.4f;
    vf3 boundary_size = { 0.6f, 0.6f, 0.6f };
    vf3 boundary_min = { -0.6f, -0.6f, -0.6f };
    vf3 boundary_max = { 0.6f, 0.6f, 0.6f };
};

namespace config
{
    // Simulation
    static const s32 NUM_PARTICLES   = 25000;
    static const f32 REST_DENSITY    = 997.0f;
    static const f32 GAS_CONSTANT    = 1.0f;     // Gas Constant, or stiffness control the particle spread, make 5 to implode
    static const f32 VISCOSITY       = 0.002f;  // Kinematic Viscosity
    static const f32 SURFACE_TENSION_CONSTANT = 0.0000001f; // 0.0728f
    static const vf3 GRAVITY         = { 0.0f, -9.80665f, 0.0f };
    static const f32 MASS            = 0.01f;
    static const f32 DT              = 0.0001f;
    // Smoothing Kernel
    static const f32 SMOOTHING_RADIUS   = 0.04f;
    static const f32 SMOOTHING_RADIUS_2 = std::pow(SMOOTHING_RADIUS, 2.0f);
    static const f32 POLY6              = 315.0f / (64.0f * PI * std::pow(SMOOTHING_RADIUS, 9.0f));   // Density
    static const f32 SPIKY_GRAD         = -45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f));           // Pressure
    static const f32 SPIKY_LAPLACIAN    =  45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f));           // Viscosity
    // Boundary
    static const f32 BOUNDARY_EPSILON = 0.0000001f;
    static const f32 BOUNDARY_DAMPING = -0.4f;
    static const vf3 BOUNDARY         = { 0.6f, 0.6f, 0.6f };
    static const vf3 BOUNDARY_MIN     = -BOUNDARY;
    static const vf3 BOUNDARY_MAX     = BOUNDARY;
    // Hash map
    const u32 TABLE_SIZE  = 262144;
    const u32 NO_PARTICLE = 0xFFFFFFFF;
}

using namespace config;

class SPH
{
public:
    struct particle
    {
        vf3 position        = { 0.0f, 0.0f, 0.0f };
        vf3 velocity        = { 0.0f, 0.0f, 0.0f };
        vf3 acceleration    = { 0.0f, 0.0f, 0.0f };

        f32 radius          = SMOOTHING_RADIUS;
        f32 mass            = MASS;
        f32 density         = REST_DENSITY;
        f32 pressure        = 0.0f;
        vf3 viscosity       = { 0.0f, 0.0f, 0.0f };
        vf3 surface_tension = { 0.0f, 0.0f, 0.0f };

        u32 hash = 0;
    };

private:
    std::vector<particle> fluid_particles;
    //std::vector<particle> boundary_particles;

    settings s;

private: // Hash map
    std::vector<std::vector<u32>> particle_table;
    vi3 cell(const particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
    u32 hash(const vi3& cell) { return ((u32)(cell.x * 73856093) ^ (u32)(cell.y * 19349663) ^ (u32)(cell.z * 83492791)) % TABLE_SIZE; }

public:
    SPH() = default;
    SPH(const settings si)
    { 
        Init(si);
    }

public:
    void Init(const settings si)
    {
        s = si;
        // Initialize Fluid Particles
        fluid_particles.clear();
        fluid_particles.reserve(s.n_particles);
        s32 grid_size      = static_cast<s32>(std::cbrt(s.n_particles));
        f32 offset         = s.smoothing_radius * 0.5f;
        f32 half_grid_size = (grid_size - 1) * offset / 2.0f;

        #pragma omp parallel
        {
            std::vector<particle> thread_local_particles;
            #pragma omp for collapse(3)
            for (s32 x = 0; x < grid_size; x++) 
            {
                for (s32 y = 0; y < grid_size; y++) 
                {
                    for (s32 z = 0; z < grid_size; z++) 
                    {
                        if (fluid_particles.size() >= s.n_particles)
                            break;

                        particle p;
                        p.position = {
                            x * offset - half_grid_size,
                            y * offset - half_grid_size + 0.25f,
                            z * offset - half_grid_size,
                        };
                        p.hash = hash(cell(p, s.smoothing_radius));

                        thread_local_particles.push_back(p);
                    }
                }
            }
            // Combine the local particles from each thread into the global vector
            #pragma omp critical
            {
                fluid_particles.insert(fluid_particles.end(), thread_local_particles.begin(), thread_local_particles.end());
            }
        }
        // Sort Particles
        std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const particle& i, const particle& j) { return i.hash < j.hash; });

        // Create Neighbor Table
        particle_table.clear();
        particle_table.resize(TABLE_SIZE);
        #pragma omp parallel
        {
            std::vector<std::vector<u32>> local_table(TABLE_SIZE);
            #pragma omp for
            for (s32 i = 0; i < fluid_particles.size(); i++) {
                u32 current_hash = fluid_particles[i].hash;
                local_table[current_hash].push_back(i);
            }
            #pragma omp critical
            {
                for (s32 i = 0; i < TABLE_SIZE; i++) {
                    particle_table[i].insert(particle_table[i].end(), local_table[i].begin(), local_table[i].end());
                }
            }
        }

        // TODO: Initialize Boundary Particles (obstacle)
    }

    void Simulate(f32 dt)
    {
        // Compute hash
        #pragma omp parallel for
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            particle& p = fluid_particles[i];
            p.hash = hash(cell(p, s.smoothing_radius));
        }

        // Sort Particles
        std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const particle& i, const particle& j) { return i.hash < j.hash; });

        // Create Neighbor Table
        particle_table.clear();
        particle_table.resize(TABLE_SIZE);

        #pragma omp parallel
        {
            std::vector<std::vector<u32>> local_table(TABLE_SIZE);
            #pragma omp for
            for (s32 i = 0; i < fluid_particles.size(); i++) {
                u32 current_hash = fluid_particles[i].hash;
                local_table[current_hash].push_back(i);
            }

            #pragma omp critical
            {
                for (s32 i = 0; i < TABLE_SIZE; i++) {
                    particle_table[i].insert(particle_table[i].end(), local_table[i].begin(), local_table[i].end());
                }
            }
        }

        //TODO: Precompute neighbor with search?

        // Dynamic Spatial Grid Search
        #pragma omp parallel for
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            particle& pi = fluid_particles[i];

            pi.density = s.rest_density;
            pi.pressure = 0.0f;

            vi3 cell_origin = cell(pi, s.smoothing_radius);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        //#pragma omp simd
                        for (s32 j = 0; j < neighbors.size(); j++)
                        {
                            const particle& pj = fluid_particles[neighbors[j]];

                            if (&pi == &pj || cell_hash != pj.hash) 
                                continue;
                            
                            f32 distance2 = glm::distance2(pj.position, pi.position);
                            if (distance2 < s.smoothing_radius2)
                            {
                                pi.density += pj.mass * s.poly6 * std::pow(s.smoothing_radius2 - distance2, 3.0f);
                            }
                        }
                    }
                }
            }
            pi.pressure = s.gas_constant * (pi.density - s.rest_density);
        }
        
        #pragma omp parallel for
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            particle& pi = fluid_particles[i];

            pi.acceleration    = {};
            vf3 pressure_force = {};
            vf3 viscous_force  = {};
            vf3 surface_tension_force = {};
            vf3 normal = {};
            vf3 gravity_force = pi.mass * s.gravity;

            vi3 cell_origin = cell(pi, s.smoothing_radius);

            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        //#pragma omp simd
                        for (s32 j = 0; j < neighbors.size(); j++)
                        {
                            const particle& pj = fluid_particles[neighbors[j]];

                            if (&pi == &pj || pj.hash != cell_hash) 
                                continue;

                            vf3 difference = pj.position - pi.position;
                            f32 distance2  = glm::length2(difference);
                            if (distance2 < s.smoothing_radius2)
                            {
                                f32 distance  = std::sqrt(distance2);
                                vf3 direction = glm::normalize(difference);

                                pressure_force += -direction * pj.mass * (pi.pressure + pj.pressure) / (2.0f * pj.density) * s.spiky_grad * std::pow(s.smoothing_radius - distance, 3.0f);
                                viscous_force  +=  pj.mass * s.viscosity * (pj.velocity - pi.velocity) / pj.density * s.spiky_laplacian * (s.smoothing_radius - distance);
                                normal         +=  direction * pj.mass / pj.density * s.spiky_grad * std::pow(s.smoothing_radius - distance, 3.0f);
                            }
                        }
                    }
                }
            }

            f32 curvature = -glm::length(normal);
            if (curvature > 0.0f)
                surface_tension_force = -s.surface_tension_constant * curvature * s.spiky_laplacian * glm::normalize(normal);

            // REVISE: Do we add the surface tension force before or after the density division?
            pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;
        }
        
        #pragma omp parallel for
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            particle& pi = fluid_particles[i];

            // Update Velocity
            pi.velocity += pi.acceleration * dt;

            // Update Position
            pi.position += pi.velocity * dt;
            
            // Particle Collision
            vi3 cell_origin = cell(pi, s.smoothing_radius);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        //#pragma omp simd
                        for (s32 j = 0; j < neighbors.size(); j++)
                        {
                            const particle& pj = fluid_particles[neighbors[j]];

                            if (&pi == &pj || pj.hash != cell_hash) 
                                continue;

                            vf3 difference = pj.position - pi.position;
                            f32 distance2  = glm::length2(difference);
                            if (distance2 < 0.0f)
                            {
                                f32 distance = std::sqrt(distance2);
                                if (distance < 1e-6f)
                                    distance = 1e-6f;

                                vf3 direction = difference / distance;
                                f32 overlap = s.smoothing_radius - distance;

                                // Move particles apart to resolve overlap
                                pi.position -= direction * overlap * 0.5f;
                                fluid_particles[j].position += direction * overlap * 0.5f;

                                // Compute relative velocity
                                vf3 relative_velocity = pi.velocity - pj.velocity;

                                // Apply impulse to separate particles
                                f32 impulse = glm::dot(relative_velocity, direction);
                                pi.velocity -= impulse * direction;
                                fluid_particles[j].velocity += impulse * direction;
                            }
                        }
                    }
                }
            }

            // Boundary Collision
            if (pi.position.x - s.boundary_epsilon < s.boundary_min.x)
            {
                pi.velocity.x *= s.boundary_damping;
                pi.position.x = s.boundary_min.x + s.boundary_epsilon;
            }
            if (pi.position.x + s.boundary_epsilon > BOUNDARY_MAX.x)
            {
                pi.velocity.x *= s.boundary_damping;
                pi.position.x = s.boundary_max.x - s.boundary_epsilon;
            }

            if (pi.position.y - s.boundary_epsilon < s.boundary_min.y)
            {
                pi.velocity.y *= s.boundary_damping;
                pi.position.y = s.boundary_min.y + s.boundary_epsilon;
            }
            if (pi.position.y + s.boundary_epsilon > s.boundary_max.y)
            {
                pi.velocity.y *= s.boundary_damping;
                pi.position.y = s.boundary_max.y - s.boundary_epsilon;
            }

            if (pi.position.z - s.boundary_epsilon < s.boundary_min.z)
            {
                pi.velocity.z *= s.boundary_damping;
                pi.position.z = s.boundary_min.z + s.boundary_epsilon;
            }
            if (pi.position.z + s.boundary_epsilon > s.boundary_max.z)
            {
                pi.velocity.z *= s.boundary_damping;
                pi.position.z = s.boundary_max.z - s.boundary_epsilon;
            }
        }
    }

public:
    std::vector<particle>& GetParticles() { return fluid_particles; }
};