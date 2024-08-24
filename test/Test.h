#pragma once

#include <map>
#include <print>
#include <chrono>
#include <algorithm>
#include <unordered_map>

#include "Common.h"
#include "Random.h"
#include "SPH.h"

// Neighbourhood search evaluation

// Particle
struct particle
{
    vf3 position = { 0.0f, 0.0f, 0.0f };
    f32 mass = 0.0f;
    f32 density = 0.0f;
    f32 pressure = 0.0f;
    u32 hash = 0;
    //std::vector<particle*> neighbors;
    std::vector<u32> neighbors;

};

// Grid-based hash neighbour search algorithm
const u32 TABLE_SIZE = 262144;
const u32 NO_PARTICLE = 0xFFFFFFFF;

static vi3 get_cell(const particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
static u32 hash(const vi3& cell) { return ((u32)(cell.x * 73856093) ^ (u32)(cell.y * 19349663) ^ (u32)(cell.z * 83492791)) % TABLE_SIZE; }

// hash_map: O(nm)
static void hash_map()
{
    std::vector<particle> particles;
    particles.reserve(config::NUM_PARTICLES);

    random rand(1337);
    for (s32 i = 0; i < config::NUM_PARTICLES; i++) {
        particle p;
        p.position = {
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
        };
        p.hash = hash(get_cell(p, config::SMOOTHING_RADIUS));  // Compute Hash
        particles.push_back(p);
    }

    // Sort Particles by Hash
    std::sort(particles.begin(), particles.end(), [&](const particle& i, const particle& j) { return i.hash < j.hash; });

    // Create Neighbor Table (storing indices of particles)
    std::vector<std::vector<u32>> particle_table(TABLE_SIZE);
    for (s32 i = 0; i < particles.size(); i++) {
        u32 current_hash = particles[i].hash;
        particle_table[current_hash].push_back(i);
    }

    // Neighbor Search
    for (auto& pi : particles) {
        pi.density = config::REST_DENSITY;
        pi.pressure = 0.0f;

        vi3 cell = get_cell(pi, config::SMOOTHING_RADIUS);
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
                        const particle& pj = particles[j];
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
}

// unordered_map: O(nm)
static void unordered_map()
{
    std::vector<particle> particles;
    particles.reserve(config::NUM_PARTICLES);

    random rand(1337);
    for (int i = 0; i < config::NUM_PARTICLES; i++) {
        particle p;
        p.position = {
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
        };
        // Compute hash
        p.hash = hash(get_cell(p, config::SMOOTHING_RADIUS));
        particles.push_back(p);
    }

    // Sort particles by their hash values
    std::sort(particles.begin(), particles.end(), [](const particle& a, const particle& b) { return a.hash < b.hash; });

    // Create neighbor table
    std::unordered_map<u32, std::vector<u32>> particle_table;
    for (u32 i = 0; i < particles.size(); i++) {
        u32 current_hash = particles[i].hash;
        particle_table[current_hash].push_back(i);
    }

    // Neighbor search
    for (auto& pi : particles) 
    {
        vi3 cell = get_cell(pi, config::SMOOTHING_RADIUS);
        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++) 
            {
                for (s32 z = -1; z <= 1; z++) 
                {
                    u32 cell_hash = hash(cell + vi3(x, y, z));
                    auto it = particle_table.find(cell_hash);
                    if (it == particle_table.end()) continue;

                    const auto& neighbor_indices = it->second;
                    for (u32 j : neighbor_indices) {
                        const particle& pj = particles[j];
                        if (&pi == &pj) continue;

                        vf3 difference = pj.position - pi.position;
                        f32 distance2 = glm::length2(difference);
                        if (distance2 < config::SMOOTHING_RADIUS_2) {
                            pi.density += pj.mass * config::POLY6 * std::pow(config::SMOOTHING_RADIUS_2 - distance2, 3.0f);
                        }
                    }
                }
            }
        }

        pi.pressure = config::GAS_CONSTANT * (pi.density - config::REST_DENSITY);
    }
}

static void precompute_neighbour()
{
    std::vector<particle> particles;
    particles.reserve(config::NUM_PARTICLES);

    random rand(1337);
    for (s32 i = 0; i < config::NUM_PARTICLES; i++)
    {
        particle p;
        p.position = {
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
        };
        // Compute Hash
        p.hash = hash(get_cell(p, config::SMOOTHING_RADIUS));
        particles.push_back(p);
    }

    std::unordered_map<u32, u32> particle_table;
    for (u32 i = 0; i < particles.size(); i++)
    {
        u32 current_hash = particles[i].hash;
        particle_table[current_hash] = i;
    }

    // Precompute neighbor
    for (auto& pi : particles)
    {
        vi3 cell = get_cell(pi, config::SMOOTHING_RADIUS);
        for (s32 x = -1; x <= 1; x++) {
            for (s32 y = -1; y <= 1; y++) {
                for (s32 z = -1; z <= 1; z++) {
                    u32 cell_hash = hash(cell + vi3(x, y, z));
                    auto it = particle_table.find(cell_hash);
                    if (it == particle_table.end()) continue;

                    const particle& pj = particles[it->second];

                    if (&pi == &pj || pj.hash != cell_hash) continue;
                    vf3 difference = pj.position - pi.position;
                    f32 distance2 = glm::length2(difference);
                    if (distance2 < config::SMOOTHING_RADIUS_2) {
                        pi.neighbors.push_back(it->second);
                    }
                }
            }
        }
    }

    // Neighbor search
    for (auto& pi : particles)
    {
        for (auto& j : pi.neighbors)
        {
            particle& pj = particles[j];
            f32 distance2 = glm::distance2(pj.position, pi.position);
            pi.density += pj.mass * config::POLY6 * std::pow(config::SMOOTHING_RADIUS_2 - distance2, 3.0f);
        }
        pi.pressure = config::GAS_CONSTANT * (pi.density - config::REST_DENSITY);
    }
}

// Two nested for loop: O(n^2)
static void naive_function()
{
    std::vector<particle> particles;

    random rand(1337);
    for (s32 i = 0; i < config::NUM_PARTICLES; i++)
    {
        particle p;
        p.position = {
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
            rand.uniform(-0.5f, 0.0f),
        };
        particles.push_back(p);
    }

    for (s32 i = 0; i < particles.size(); i++)
    {
        particle& pi = particles[i];
        pi.density = config::REST_DENSITY;
        pi.pressure = 0.0f;

        for (s32 j = 0; j < particles.size(); j++)
        {
            particle& pj = particles[j];

            if (&pi == &pj) continue;

            f32 distance2 = glm::distance2(pj.position, pi.position);
            if (distance2 <= config::SMOOTHING_RADIUS_2)
            {
                pi.density += pj.mass * config::POLY6 * std::pow(config::SMOOTHING_RADIUS_2 - distance2, 3.0f);
            }
        }
        pi.pressure = config::GAS_CONSTANT * (pi.density - config::REST_DENSITY);
    }
}

static void naive_simulate()
{
/*
    // Neighbourhood Search
    // TODO: Decouple with density and pressute calculation
    // Naive search implementation
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        particle& pi = fluid_particles[i];
        pi.density = config::REST_DENSITY;
        pi.pressure = 0.0f;

        for (s32 j = 0; j < fluid_particles.size(); j++)
        {
            particle& pj = fluid_particles[j];

            if (&pi == &pj) continue;

            f32 distance2 = glm::distance2(pj.position, pi.position);
            if (distance2 <= config::SMOOTHING_RADIUS_2)
            {
                pi.density += pj.mass * config::POLY6 * std::pow(config::SMOOTHING_RADIUS_2 - distance2, 3.0f);
            }
        }
        pi.pressure = config::GAS_CONSTANT * (pi.density - config::REST_DENSITY);
    }

    // Acceleration
        // Internal Force = Pressure Force + Viscous Force
        // External Force = Surface Tension Force + Gravity Force
        // Acceleration   = Internal Force + External Force
    // TODO: Implement Versatile Surface Tension and Adhesion in SPH Fluids

    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        particle& pi = fluid_particles[i];
        pi.acceleration           = {};
        vf3 pressure_force        = {};
        vf3 viscous_force         = {};
        vf3 surface_tension_force = {};
        vf3 normal                = {};
        vf3 gravity_force         = pi.mass * config::GRAVITY;

        for (s32 j = 0; j < fluid_particles.size(); j++)
        {
            particle& pj = fluid_particles[j];

            if (&pi == &pj) continue;

            vf3 difference = pj.position - pi.position;
            f32 distance2  = glm::length2(difference);
            if (distance2 < config::SMOOTHING_RADIUS_2)
            {
                f32 distance  = std::sqrt(distance2);
                vf3 direction = glm::normalize(difference);

                pressure_force += -direction * pj.mass * (pi.pressure + pj.pressure) / (2.0f * pj.density) * config::SPIKY * std::pow(config::SMOOTHING_RADIUS - distance, 3.0f);
                viscous_force  += pj.mass * config::VISCOSITY * (pj.velocity - pi.velocity) / pj.density * config::LAPLACIAN * (config::SMOOTHING_RADIUS - distance);

                normal         += direction * pj.mass / pj.density * config::SPIKY * std::pow(config::SMOOTHING_RADIUS - distance, 2.0f);
            }
        }

        f32 curvature = -glm::length(normal);
        if (curvature > 0.0f)
            surface_tension_force = -config::SURFACE_TENSION_CONSTANT * curvature * config::LAPLACIAN * glm::normalize(normal);

        // REVISE: Do we add the surface tension force before or after the density division?
        pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;
    }

    // Integrate
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        particle& pi = fluid_particles[i];

        // Update Velocity
        pi.velocity += pi.acceleration * dt;

        // Update Position
        pi.position += pi.velocity * dt;

        // Particle Collsion
        for (s32 j = 0; j < fluid_particles.size(); j++)
        {
            particle& pj = fluid_particles[j];

            if (&pi == &pj) continue;

            vf3 difference = pj.position - pi.position;
            f32 distance2  = glm::length2(difference);
            if (distance2 < 0.0f)
            {
                f32 distance = std::sqrt(distance2);
                vf3 surface_normal = -difference / distance;
                pi.position = pi.position + (difference / distance) * config::SMOOTHING_RADIUS;
                pi.velocity = pi.velocity - difference * 2.0f * glm::dot(pi.velocity, surface_normal);
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
*/
}