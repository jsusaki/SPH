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
namespace test {
    // Particle
    struct Particle
    {
        vf3 position = { 0.0f, 0.0f, 0.0f };
        f32 mass = 0.0f;
        f32 density = 0.0f;
        f32 pressure = 0.0f;

        u32 hash = 0;
        std::vector<u32> neighbors;
    };

    // Grid-based hash neighbor search algorithm
    const u32 TABLE_SIZE = 262144;
    const u32 NO_PARTICLE = 0xFFFFFFFF;
    static vi3 cell(const Particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
    static u32 hash(const vi3& cell) { return ((u32)(cell.x * 73856093) ^ (u32)(cell.y * 19349663) ^ (u32)(cell.z * 83492791)) % TABLE_SIZE; }

    // hash_map: O(nm)
    static void spatial_hash_map()
    {
        std::vector<Particle> particles;
        particles.reserve(NUM_PARTICLES);

        random rand(1337);
        for (s32 i = 0; i < NUM_PARTICLES; i++) {
            Particle p;
            p.position = {
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
            };
            p.hash = hash(cell(p, SMOOTHING_RADIUS));  // Compute Hash
            particles.push_back(p);
        }

        // Sort Particles by Hash
        std::sort(particles.begin(), particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

        // Create Neighbor Table (storing indices of particles)
        std::vector<std::vector<u32>> particle_table(TABLE_SIZE);
        for (s32 i = 0; i < particles.size(); i++) {
            u32 current_hash = particles[i].hash;
            particle_table[current_hash].push_back(i);
        }

        // Neighborhood Search
        for (auto& pi : particles) 
        {
            pi.density = REST_DENSITY;
            pi.pressure = 0.0f;

            vi3 cell_origin = cell(pi, SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        for (u32 j : neighbors)
                        {
                            const Particle& pj = particles[j];
                            if (&pi == &pj || cell_hash != pj.hash) continue;
                            f32 distance2 = glm::distance2(pj.position, pi.position);
                            if (distance2 < SMOOTHING_RADIUS_2)
                            {
                                pi.density += pj.mass * POLY6 * std::pow(SMOOTHING_RADIUS_2 - distance2, 3.0f);
                            }
                        }
                    }
                }
            }
            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
        }
    }

    static void intrinsic_spatial_hash_map()
    {
        std::vector<Particle> particles;
        particles.reserve(NUM_PARTICLES);

        random rand(1337);
        for (s32 i = 0; i < NUM_PARTICLES; i++) {
            Particle p;
            p.position = {
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
            };
            p.hash = hash(cell(p, SMOOTHING_RADIUS));  // Compute Hash
            particles.push_back(p);
        }

        // Sort Particles by Hash
        std::sort(particles.begin(), particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

        // Create Neighbor Table (storing indices of particles)
        std::vector<std::vector<u32>> particle_table(TABLE_SIZE);
        for (s32 i = 0; i < particles.size(); i++) 
        {
            u32 current_hash = particles[i].hash;
            particle_table[current_hash].push_back(i);
        }

        f32 pjx_array[8] = { 0.0f }, pjy_array[8] = { 0.0f }, pjz_array[8] = { 0.0f };

        // 64-bit float registers
        // Positions
        __m256 _pjx, _pjy, _pjz;
        __m256 _pix, _piy, _piz;

        // Squared distance
        __m256 _dx, _dy, _dz, _dx2, _dy2, _dz2, _dx2dy2;
        __m256 _dist2;
        // if condition with bitwise AND operation
        __m256 _mask;
        // Density with poly6
        __m256 _poly6, _smoothing_radius2, _mass;
        __m256 _radius_diff, _radius_diff2, _radius_diff3, _poly6_radius_diff;
        __m256 _density;
        __m128 _lo, _hi;
        __m128 _sum;

        // Set constants into AVX2 registers
        _poly6 = _mm256_set1_ps(POLY6);
        _smoothing_radius2 = _mm256_set1_ps(SMOOTHING_RADIUS_2);
        _mass = _mm256_set1_ps(MASS);

        for (s32 i = 0; i < particles.size(); i++)
        {
            Particle& pi = particles[i];
            pi.density = REST_DENSITY;
            pi.pressure = 0.0f;

            // Neighborhood Search
            vi3 cell_origin = cell(pi, SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];

                        for (s32 j = 0; j < neighbors.size(); j += 8) // Process 8 neighbors at a time
                        {
                            // Transform into a contiguous array for SIMD
                            for (s32 k = 0; k < 8; k++) 
                            {
                                if (j + k < neighbors.size()) 
                                {
                                    pjx_array[k] = particles[neighbors[j + k]].position.x;
                                    pjy_array[k] = particles[neighbors[j + k]].position.y;
                                    pjz_array[k] = particles[neighbors[j + k]].position.z;
                                }
                            }

                            // Load neighbor positions
                            _pjx = _mm256_load_ps(pjx_array);
                            _pjy = _mm256_load_ps(pjy_array);
                            _pjz = _mm256_load_ps(pjz_array);
                            // Set own positions
                            _pix = _mm256_set1_ps(pi.position.x);
                            _piy = _mm256_set1_ps(pi.position.y);
                            _piz = _mm256_set1_ps(pi.position.z);

                            // Compute squared distance
                            // f32 distance2 = glm::distance2(pj.position, pi.position);
                            _dx     = _mm256_sub_ps(_pix, _pjx);
                            _dy     = _mm256_sub_ps(_piy, _pjy);
                            _dz     = _mm256_sub_ps(_piz, _pjz);
                            _dx2    = _mm256_mul_ps(_dx, _dx);
                            _dy2    = _mm256_mul_ps(_dy, _dy);
                            _dz2    = _mm256_mul_ps(_dz, _dz);
                            _dx2dy2 = _mm256_add_ps(_dx2, _dy2);
                            _dist2  = _mm256_add_ps(_dx2dy2, _dz2);

                            // Mask particles within squared smoothing radius
                            // if (distance2 < settings.smoothing_radius2)
                            _mask = _mm256_cmp_ps(_dist2, _smoothing_radius2, _CMP_LT_OS);

                            // Compute density for particles within the smoothing radius
                            //pi.density += pj.mass * settings.poly6 * std::pow(settings.smoothing_radius2 - distance2, 3.0f);
                            _radius_diff = _mm256_sub_ps(_smoothing_radius2, _dist2);
                            _radius_diff2 = _mm256_mul_ps(_radius_diff, _radius_diff);
                            _radius_diff3 = _mm256_mul_ps(_radius_diff2, _radius_diff);
                            _poly6_radius_diff = _mm256_mul_ps(_poly6, _radius_diff3);
                            _density = _mm256_mul_ps(_mass, _poly6_radius_diff);
                            // Apply mask to zero out the contribution for particles outside the smoothing radius
                            _density = _mm256_blendv_ps(_mm256_setzero_ps(), _density, _mask);

                            // Reduce
                            // Horizontally add pairs of adjacent elements
                            _lo = _mm256_castps256_ps128(_density);
                            _hi = _mm256_extractf128_ps(_density, 1);
                            _lo = _mm_add_ps(_lo, _hi);
                            // Horizontally add the resulting vector
                            _sum = _mm_hadd_ps(_lo, _lo);
                            _sum = _mm_hadd_ps(_sum, _sum);
                            // Horizontally add to sum the contributions
                            pi.density += _mm_cvtss_f32(_sum);
                        }
                    }
                }
            }
            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
        }
    }

    static void parallel_spatial_hash_map()
    {
        std::vector<Particle> particles;
        particles.reserve(NUM_PARTICLES);
        random rand(1337);
        for (s32 i = 0; i < NUM_PARTICLES; i++) 
        {
            Particle p;
            p.position = {
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
            };
            p.hash = hash(cell(p, SMOOTHING_RADIUS));  // Compute Hash
            particles.push_back(p);
        }

        // Compute hash
        #pragma omp parallel for
        for (s32 i = 0; i < particles.size(); i++)
        {
            Particle& p = particles[i];
            p.hash = hash(cell(p, SMOOTHING_RADIUS));
        }

        // Sort Particles
        std::sort(particles.begin(), particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

        std::vector<std::vector<u32>> particle_table;
        // Create Neighbor Table
        particle_table.clear();
        particle_table.resize(TABLE_SIZE);

        #pragma omp parallel
        {
            std::vector<std::vector<u32>> local_table(TABLE_SIZE);
            #pragma omp for
            for (s32 i = 0; i < particles.size(); i++) 
            {
                u32 current_hash = particles[i].hash;
                local_table[current_hash].push_back(i);
            }

            #pragma omp critical
            {
                for (s32 i = 0; i < TABLE_SIZE; i++) 
                {
                    particle_table[i].insert(particle_table[i].end(), local_table[i].begin(), local_table[i].end());
                }
            }
        }

        // Dynamic Spatial Grid Search
        #pragma omp parallel for
        for (s32 i = 0; i < particles.size(); i++)
        {
            Particle& pi = particles[i];

            pi.density = REST_DENSITY;
            pi.pressure = 0.0f;

            vi3 cell_origin = cell(pi, SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        const auto& neighbors = particle_table[cell_hash];
                        for (u32 j : neighbors)
                        {
                            const Particle& pj = particles[j];

                            if (&pi == &pj || cell_hash != pj.hash)
                                continue;

                            f32 distance2 = glm::distance2(pj.position, pi.position);
                            if (distance2 < SMOOTHING_RADIUS_2)
                            {
                                pi.density += pj.mass * POLY6 * std::pow(SMOOTHING_RADIUS_2 - distance2, 3.0f);
                            }
                        }
                    }
                }
            }
            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
        }
    }

    // unordered_map: O(nm)
    static void unordered_map()
    {
        std::vector<Particle> particles;
        particles.reserve(NUM_PARTICLES);

        random rand(1337);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            Particle p;
            p.position = {
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
            };
            // Compute hash
            p.hash = hash(cell(p, SMOOTHING_RADIUS));
            particles.push_back(p);
        }

        // Sort particles by their hash values
        std::sort(particles.begin(), particles.end(), [](const Particle& a, const Particle& b) { return a.hash < b.hash; });

        // Create neighbor table
        std::unordered_map<u32, std::vector<u32>> particle_table;
        for (u32 i = 0; i < particles.size(); i++) {
            u32 current_hash = particles[i].hash;
            particle_table[current_hash].push_back(i);
        }

        // Neighbor search
        for (auto& pi : particles)
        {
            vi3 cell_origin = cell(pi, SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++)
            {
                for (s32 y = -1; y <= 1; y++)
                {
                    for (s32 z = -1; z <= 1; z++)
                    {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        auto it = particle_table.find(cell_hash);
                        if (it == particle_table.end()) continue;

                        for (u32 j : it->second) {
                            const Particle& pj = particles[j];
                            if (&pi == &pj || pj.hash != cell_hash) continue;

                            vf3 difference = pj.position - pi.position;
                            f32 distance2 = glm::length2(difference);
                            if (distance2 < SMOOTHING_RADIUS_2) {
                                pi.density += pj.mass * POLY6 * std::pow(SMOOTHING_RADIUS_2 - distance2, 3.0f);
                            }
                        }
                    }
                }
            }

            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
        }
    }

    static void precompute_neighbor()
    {
        std::vector<Particle> particles;
        particles.reserve(NUM_PARTICLES);

        random rand(1337);
        for (s32 i = 0; i < NUM_PARTICLES; i++)
        {
            Particle p;
            p.position = {
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
            };
            // Compute Hash
            p.hash = hash(cell(p, SMOOTHING_RADIUS));
            particles.push_back(p);
        }

        // Create Neighbor Table (storing indices of particles)
        std::unordered_map<u32, std::vector<u32>> particle_table;
        for (u32 i = 0; i < particles.size(); i++) {
            u32 current_hash = particles[i].hash;
            particle_table[current_hash].push_back(i);
        }

        // Precompute neighbor
        for (auto& pi : particles)
        {
            vi3 cell_origin = cell(pi, SMOOTHING_RADIUS);
            for (s32 x = -1; x <= 1; x++) {
                for (s32 y = -1; y <= 1; y++) {
                    for (s32 z = -1; z <= 1; z++) {
                        u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                        auto it = particle_table.find(cell_hash);
                        if (it == particle_table.end()) continue;

                        //(something is sus here...)
                        for (u32 j : it->second) {
                            const Particle& pj = particles[j];
                            if (&pi == &pj || pj.hash != cell_hash) continue;
                            vf3 difference = pj.position - pi.position;
                            f32 distance2 = glm::length2(difference);
                            if (distance2 < SMOOTHING_RADIUS_2) {
                                pi.neighbors.push_back(j);
                            }
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
                Particle& pj = particles[j];
                f32 distance2 = glm::distance2(pj.position, pi.position);
                pi.density += pj.mass * POLY6 * std::pow(SMOOTHING_RADIUS_2 - distance2, 3.0f);
            }
            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
        }
    }

    // Two nested for loop: O(n^2)
    static void naive()
    {
        std::vector<Particle> particles;

        random rand(1337);
        for (s32 i = 0; i < NUM_PARTICLES; i++)
        {
            Particle p;
            p.position = {
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
                rand.uniform(-0.5f, 0.0f),
            };
            particles.push_back(p);
        }

        for (s32 i = 0; i < particles.size(); i++)
        {
            Particle& pi = particles[i];
            pi.density = REST_DENSITY;
            pi.pressure = 0.0f;

            for (s32 j = 0; j < particles.size(); j++)
            {
                Particle& pj = particles[j];

                if (&pi == &pj) continue;

                f32 distance2 = glm::distance2(pj.position, pi.position);
                if (distance2 <= SMOOTHING_RADIUS_2)
                {
                    pi.density += pj.mass * POLY6 * std::pow(SMOOTHING_RADIUS_2 - distance2, 3.0f);
                }
            }
            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
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
            pi.density = REST_DENSITY;
            pi.pressure = 0.0f;

            for (s32 j = 0; j < fluid_particles.size(); j++)
            {
                particle& pj = fluid_particles[j];

                if (&pi == &pj) continue;

                f32 distance2 = glm::distance2(pj.position, pi.position);
                if (distance2 <= SMOOTHING_RADIUS_2)
                {
                    pi.density += pj.mass * POLY6 * std::pow(SMOOTHING_RADIUS_2 - distance2, 3.0f);
                }
            }
            pi.pressure = GAS_CONSTANT * (pi.density - REST_DENSITY);
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
            vf3 gravity_force         = pi.mass * GRAVITY;

            for (s32 j = 0; j < fluid_particles.size(); j++)
            {
                particle& pj = fluid_particles[j];

                if (&pi == &pj) continue;

                vf3 difference = pj.position - pi.position;
                f32 distance2  = glm::length2(difference);
                if (distance2 < SMOOTHING_RADIUS_2)
                {
                    f32 distance  = std::sqrt(distance2);
                    vf3 direction = glm::normalize(difference);

                    pressure_force += -direction * pj.mass * (pi.pressure + pj.pressure) / (2.0f * pj.density) * SPIKY * std::pow(SMOOTHING_RADIUS - distance, 3.0f);
                    viscous_force  += pj.mass * VISCOSITY * (pj.velocity - pi.velocity) / pj.density * LAPLACIAN * (SMOOTHING_RADIUS - distance);

                    normal         += direction * pj.mass / pj.density * SPIKY * std::pow(SMOOTHING_RADIUS - distance, 2.0f);
                }
            }

            f32 curvature = -glm::length(normal);
            if (curvature > 0.0f)
                surface_tension_force = -SURFACE_TENSION_CONSTANT * curvature * LAPLACIAN * glm::normalize(normal);

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
                    pi.position = pi.position + (difference / distance) * SMOOTHING_RADIUS;
                    pi.velocity = pi.velocity - difference * 2.0f * glm::dot(pi.velocity, surface_normal);
                }
            }

            // Boundary Collision
            if (pi.position.x - BOUNDARY_EPS < BOUNDARY_MIN.x)
            {
                pi.velocity.x *= BOUNDARY_DAMPING;
                pi.position.x = BOUNDARY_MIN.x + BOUNDARY_EPS;
            }
            if (pi.position.x + BOUNDARY_EPS > BOUNDARY_MAX.x)
            {
                pi.velocity.x *= BOUNDARY_DAMPING;
                pi.position.x = BOUNDARY_MAX.x - BOUNDARY_EPS;
            }

            if (pi.position.y - BOUNDARY_EPS < BOUNDARY_MIN.y)
            {
                pi.velocity.y *= BOUNDARY_DAMPING;
                pi.position.y = BOUNDARY_MIN.y + BOUNDARY_EPS;
            }
            if (pi.position.y + BOUNDARY_EPS > BOUNDARY_MAX.y)
            {
                pi.velocity.y *= BOUNDARY_DAMPING;
                pi.position.y = BOUNDARY_MAX.y - BOUNDARY_EPS;
            }

            if (pi.position.z - BOUNDARY_EPS < BOUNDARY_MIN.z)
            {
                pi.velocity.z *= BOUNDARY_DAMPING;
                pi.position.z = BOUNDARY_MIN.z + BOUNDARY_EPS;
            }
            if (pi.position.z + BOUNDARY_EPS > BOUNDARY_MAX.z)
            {
                pi.velocity.z *= BOUNDARY_DAMPING;
                pi.position.z = BOUNDARY_MAX.z - BOUNDARY_EPS;
            }
        }
    */
    }


    static void openmp()
    {
        /*
        // OpenMP minimal example with multi-threading and SIMD
        #include <cmath>
        #include <omp.h>
        #include <iostream>
        int main()
        {
            const int size = 256;
            double sinTable[size];

            #pragma omp parallel for simd
            for (int n = 0; n < size; ++n)
            {
                sinTable[n] = std::sin(2 * PI * n / size);
                std::cout << "Thread " << omp_get_thread_num() << " is working on iteration " << n << "\n";
            }
        }
        // compile: g++ main.cpp -fopenmp
        */
    }
}