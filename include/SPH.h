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
            Pressure : Spikey Kernel
            Viscosity: Laplacian Kernel

        Density Function
        Pressure Function
        Viscosity Function
        Surface Tension Function
        Gravity Function
        Acceleration Function
        Velocity Function
        Position Function

        Neighbourhood Search
            Grid-based 
            Hash Table

        Time Integration (Euler Integration)

    ----------
    References
    ----------

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

*/
#pragma once

#include <vector>
#include <print>

#include "Common.h"
#include "Random.h"

const u32 SCREEN_WIDTH  = 800;
const u32 SCREEN_HEIGHT = 600;

static const f32 PI = 3.14159265358979323846f;

// TODO: make struct
namespace config
{
    // Simulation
    static const s32 NUM_PARTICLES   = 5000;
    static const f32 REST_DENSITY    = 997.0f;
    static const f32 GAS_CONSTANT    = 1.0f;     // Gas Constant, or stiffness control the particle spread, make 5 to implode
    static const f32 VISCOSITY       = 0.0005f;
    static const f32 SURFACE_TENSION_CONSTANT = 0.01f; // 0.0728f
    static const vf3 GRAVITY         = { 0.0f, -9.80665f, 0.0f };
    static const f32 MASS            = 0.01f;
    static const f32 DT              = 0.0001f;

    // Smoothing Kernel
    static const f32 SMOOTHING_RADIUS   = 0.04f;
    static const f32 SMOOTHING_RADIUS_2 = std::pow(SMOOTHING_RADIUS, 2.0f);
    //static const f32 SUPPORT_RADIUS     = 0.1f;
    static const f32 POLY6              = 315.0f / (64.0f * PI * std::pow(SMOOTHING_RADIUS, 9.0f));   // Density
    static const f32 SPIKY              = -45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f));           // Pressure
    static const f32 LAPLACIAN          =  45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f));           // Viscosity

    // Boundary
    static const f32 BOUNDARY_EPS     = 0.0000001f;
    static const f32 BOUNDARY_DAMPING = -0.3f;
    static const f32 BOUNDARY_WIDTH   = 0.5f;
    static const f32 BOUNDARY_HEIGHT  = 0.5f;
    static const f32 BOUNDARY_DEPTH   = 0.5f;
    static const vf3 BOUNDARY_MIN     = {-BOUNDARY_WIDTH,-BOUNDARY_HEIGHT,-BOUNDARY_DEPTH };
    static const vf3 BOUNDARY_MAX     = { BOUNDARY_WIDTH, BOUNDARY_HEIGHT, BOUNDARY_DEPTH };

    // Grid-based hash neighbour search algorithm
    const u32 TABLE_SIZE = 262144;
    const u32 NO_PARTICLE = 0xFFFFFFFF;
}


class SPH
{
public:
    // Particle
    struct particle
    {
        vf3 position        = { 0.0f, 0.0f, 0.0f };
        vf3 velocity        = { 0.0f, 0.0f, 0.0f };
        vf3 acceleration    = { 0.0f, 0.0f, 0.0f };

        f32 radius          = config::SMOOTHING_RADIUS;
        f32 mass            = config::MASS;
        f32 density         = config::REST_DENSITY;
        f32 pressure        = 0.0f;
        vf3 viscosity       = { 0.0f, 0.0f, 0.0f };
        vf3 surface_tension = { 0.0f, 0.0f, 0.0f };

        u32 hash = 0;
    };

private: // Hash map
    std::vector<std::vector<u32>> particle_table;
    vi3 cell(const particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
    u32 hash(const vi3& cell) { return ((u32)(cell.x * 73856093) ^ (u32)(cell.y * 19349663) ^ (u32)(cell.z * 83492791)) % config::TABLE_SIZE; }

private:
    std::vector<particle> fluid_particles;
    std::vector<particle> boundary_particles;

public:
    SPH() { Init(); }

public:
    void Init()
    {
        // Initialize Fluid Particles
        fluid_particles.clear();

        s32 grid_size = static_cast<s32>(std::cbrt(config::NUM_PARTICLES));
        f32 offset = 0.03f;
        f32 half_grid_size = (grid_size - 1) * offset / 2.0f;

        for (s32 x = 0; x < grid_size; x++)
            for (s32 y = 0; y < grid_size; y++)
                for (s32 z = 0; z < grid_size; z++)
                {
                    if (fluid_particles.size() >= config::NUM_PARTICLES) 
                        break;

                    particle p;
                    p.position = {
                        x * offset - half_grid_size,
                        y * offset - half_grid_size + 0.25f,
                        z * offset - half_grid_size,
                    };
                    p.hash = hash(cell(p, config::SMOOTHING_RADIUS));
                    fluid_particles.push_back(p);
                }

        // Sort Particles
        std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const particle& i, const particle& j) { return i.hash < j.hash; });

        // Create Neighbor Table
        particle_table = std::vector<std::vector<u32>>(config::TABLE_SIZE);
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            u32 current_hash = fluid_particles[i].hash;
            particle_table[current_hash].push_back(i);
        }

        // TODO: Initialize Boundary Particles (obstacle)
    }

    void RandomInit()
    {
        random rand(1337);
        for (s32 i = 0; i < config::NUM_PARTICLES; i++)
        {
            particle p;
            p.position = {
                rand.uniform(0.0f, 0.5f),
                rand.uniform(0.0f, 0.5f),
                rand.uniform(0.0f, 0.5f),
            };

            fluid_particles.push_back(p);
        }
    }

    void Simulate(f32 dt)
    {
        // Compute hash
        for (auto& p : fluid_particles)
            p.hash = hash(cell(p, config::SMOOTHING_RADIUS));

        // Sort Particles
        std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const particle& i, const particle& j) { return i.hash < j.hash; });

        // Create Neighbor Table
        particle_table = std::vector<std::vector<u32>>(config::TABLE_SIZE);
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            u32 current_hash = fluid_particles[i].hash;
            particle_table[current_hash].push_back(i);
        }

        // Dynamic Grid Search
        for (auto& pi : fluid_particles)
        {
            pi.density = config::REST_DENSITY;
            pi.pressure = 0.0f;

            vi3 cell = cell(pi, config::SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        for (u32 j : neighbors)
                        {
                            const particle& pj = fluid_particles[j];
                            if (&pi == &pj || cell_hash != pj.hash) continue;
                            f32 distance2 = glm::distance2(pj.position, pi.position);
                            if (distance2 < config::SMOOTHING_RADIUS_2)
                            {
                                pi.density += pj.mass * config::POLY6 * std::pow(config::SMOOTHING_RADIUS_2 - distance2, 3.0f);
                            }
                        }
                    }
                }
            }
            pi.pressure = config::GAS_CONSTANT * (pi.density - config::REST_DENSITY);
        }

        for (auto& pi : fluid_particles)
        {
            pi.acceleration = {};
            vf3 pressure_force = {};
            vf3 viscous_force = {};
            vf3 surface_tension_force = {};
            vf3 normal = {};
            vf3 gravity_force = pi.mass * config::GRAVITY;

            vi3 cell = cell(pi, config::SMOOTHING_RADIUS);

            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        for (u32 j : neighbors)
                        {
                            const particle& pj = fluid_particles[j];
                            if (&pi == &pj || pj.hash != cell_hash) continue;

                            vf3 difference = pj.position - pi.position;
                            f32 distance2 = glm::length2(difference);
                            if (distance2 < config::SMOOTHING_RADIUS_2)
                            {
                                f32 distance = std::sqrt(distance2);
                                vf3 direction = glm::normalize(difference);

                                pressure_force += -direction * pj.mass * (pi.pressure + pj.pressure) / (2.0f * pj.density) * config::SPIKY * std::pow(config::SMOOTHING_RADIUS - distance, 3.0f);
                                viscous_force += pj.mass * config::VISCOSITY * (pj.velocity - pi.velocity) / pj.density * config::LAPLACIAN * (config::SMOOTHING_RADIUS - distance);

                                normal += direction * pj.mass / pj.density * config::SPIKY * std::pow(config::SMOOTHING_RADIUS - distance, 2.0f);
                            }
                        }

                    }
                }
            }

            f32 curvature = -glm::length(normal);
            if (curvature > 0.0f)
                surface_tension_force = -config::SURFACE_TENSION_CONSTANT * curvature * config::LAPLACIAN * glm::normalize(normal);

            // REVISE: Do we add the surface tension force before or after the density division?
            pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;
        }

        for (auto& pi : fluid_particles)
        {
            // Update Velocity
            pi.velocity += pi.acceleration * dt;

            // Update Position
            pi.position += pi.velocity * dt;
            
            // Particle Collision
            vi3 cell = cell(pi, config::SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        for (u32 j : neighbors)
                        {
                            const particle& pj = fluid_particles[j];
                            if (&pi == &pj || pj.hash != cell_hash) continue;
                            vf3 difference = pj.position - pi.position;
                            f32 distance2 = glm::length2(difference);
                            if (distance2 < 0.0f)
                            {
                                f32 distance = std::sqrt(distance2);
                                vf3 surface_normal = -difference / distance;
                                pi.position = pi.position + (difference / distance) * config::SMOOTHING_RADIUS;
                                pi.velocity = pi.velocity - difference * 2.0f * glm::dot(pi.velocity, surface_normal);
                            }
                        }
                    }
                }
            }

            // Boundary Collision
            if (pi.position.x - config::BOUNDARY_EPS < config::BOUNDARY_MIN.x)
            {
                pi.velocity.x *= config::BOUNDARY_DAMPING;
                pi.position.x = config::BOUNDARY_MIN.x + config::BOUNDARY_EPS;
            }
            if (pi.position.x + config::BOUNDARY_EPS > config::BOUNDARY_MAX.x)
            {
                pi.velocity.x *= config::BOUNDARY_DAMPING;
                pi.position.x = config::BOUNDARY_MAX.x - config::BOUNDARY_EPS;
            }

            if (pi.position.y - config::BOUNDARY_EPS < config::BOUNDARY_MIN.y)
            {
                pi.velocity.y *= config::BOUNDARY_DAMPING;
                pi.position.y = config::BOUNDARY_MIN.y + config::BOUNDARY_EPS;
            }
            if (pi.position.y + config::BOUNDARY_EPS > config::BOUNDARY_MAX.y)
            {
                pi.velocity.y *= config::BOUNDARY_DAMPING;
                pi.position.y = config::BOUNDARY_MAX.y - config::BOUNDARY_EPS;
            }

            if (pi.position.z - config::BOUNDARY_EPS < config::BOUNDARY_MIN.z)
            {
                pi.velocity.z *= config::BOUNDARY_DAMPING;
                pi.position.z = config::BOUNDARY_MIN.z + config::BOUNDARY_EPS;
            }
            if (pi.position.z + config::BOUNDARY_EPS > config::BOUNDARY_MAX.z)
            {
                pi.velocity.z *= config::BOUNDARY_DAMPING;
                pi.position.z = config::BOUNDARY_MAX.z - config::BOUNDARY_EPS;
            }
        }
    }

public:
    std::vector<particle>& GetParticles() { return fluid_particles; }
};