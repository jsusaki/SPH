#include "SPH.h"


SPH::SPH() {}

SPH::SPH(const SPHSettings& s) { Init(s); }

void SPH::Init(const SPHSettings& s)
{
    settings = s;
    // Initialize Fluid Particles
    fluid_particles.clear();
    fluid_particles.reserve(settings.n_particles);
    s32 grid_size = static_cast<s32>(std::cbrt(settings.n_particles));
    std::printf("grid_size=%d\n", grid_size);
    f32 offset = settings.smoothing_radius * 0.5f;
    f32 half_grid_size = (grid_size - 1) * offset / 2.0f;

    for (s32 x = 0; x < grid_size; x++)
    {
        for (s32 y = 0; y < grid_size; y++)
        {
            for (s32 z = 0; z < grid_size; z++)
            {
                if (fluid_particles.size() >= settings.n_particles)
                    break;

                Particle p;
                p.position = {
                    x * offset - half_grid_size,
                    y * offset - half_grid_size + 0.25f,
                    z * offset - half_grid_size,
                };
                p.velocity     = { 0.0f, 0.0f, 0.0f };
                p.acceleration = { 0.0f, 0.0f, 0.0f };
                p.mass         = s.mass;
                p.pressure     = 0.0f;
                p.density      = s.rest_density;
                p.radius       = s.smoothing_radius;
                p.viscosity    = s.viscosity;
                p.hash         = hash(cell(p, settings.smoothing_radius));

                fluid_particles.push_back(p);
            }
        }
    }

    // Sort Particles
    std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

    CreateNeighborTable();

    // TODO: Initialize Boundary Particles (obstacle)

}


void SPH::Simulate(f32 dt)
{
    // Compute hash
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& p = fluid_particles[i];
        p.hash = hash(cell(p, settings.smoothing_radius));
    }
    // Sort Particles
    std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

    // Create Neighbor Table
    CreateNeighborTable();

    //TODO: Precompute neighbor with search?

    //ComputeDensityPressure();
    ComputeDensityPressureIntrinsics();

    ComputePressureViscousSurfaceTensionForce();
    //ComputePressureViscousSurfaceTensionForceIntrinsics();

    Integrate(dt);

    ComputeBoundaryCondition();
}

void SPH::ComputeDensityPressure()
{
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        pi.density = settings.rest_density;
        pi.pressure = 0.0f;

        // Spatial Hash Grid Search
        vi3 cell_origin = cell(pi, settings.smoothing_radius);
        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++)
            {
                for (s32 z = -1; z <= 1; z++)
                {
                    u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                    const auto& neighbors = particle_table[cell_hash];
                    for (s32 j = 0; j < neighbors.size(); j++)
                    {
                        Particle& pj = fluid_particles[neighbors[j]];

                        if (&pi == &pj || cell_hash != pj.hash)
                            continue;

                        f32 distance2 = glm::distance2(pj.position, pi.position);
                        if (distance2 < settings.smoothing_radius2)
                        {
                            pi.density += pj.mass * settings.poly6 * std::pow(settings.smoothing_radius2 - distance2, 3.0f);
                        }
                    }
                }
            }
        }
        pi.pressure = settings.gas_constant * (pi.density - settings.rest_density);
    }
}

void SPH::ComputeDensityPressureIntrinsics()
{
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
    _poly6 = _mm256_set1_ps(settings.poly6);
    _smoothing_radius2 = _mm256_set1_ps(settings.smoothing_radius2);
    _mass = _mm256_set1_ps(settings.mass);

    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];
        pi.density = settings.rest_density;
        pi.pressure = 0.0f;

        // Set own positions
        _pix = _mm256_set1_ps(pi.position.x);
        _piy = _mm256_set1_ps(pi.position.y);
        _piz = _mm256_set1_ps(pi.position.z);

        vi3 cell_origin = cell(pi, settings.smoothing_radius);
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
                        for (s32 k = 0; k < 8; k++) {
                            if (j + k < neighbors.size()) {
                                pjx_array[k] = fluid_particles[neighbors[j + k]].position.x;
                                pjy_array[k] = fluid_particles[neighbors[j + k]].position.y;
                                pjz_array[k] = fluid_particles[neighbors[j + k]].position.z;
                            }
                        }

                        // Load neighbor positions
                        _pjx = _mm256_load_ps(pjx_array);
                        _pjy = _mm256_load_ps(pjy_array);
                        _pjz = _mm256_load_ps(pjz_array);

                        // Compute squared distance
                        _dx     = _mm256_sub_ps(_pix, _pjx);
                        _dy     = _mm256_sub_ps(_piy, _pjy);
                        _dz     = _mm256_sub_ps(_piz, _pjz);
                        _dx2    = _mm256_mul_ps(_dx, _dx);
                        _dy2    = _mm256_mul_ps(_dy, _dy);
                        _dz2    = _mm256_mul_ps(_dz, _dz);
                        _dx2dy2 = _mm256_add_ps(_dx2, _dy2);
                        _dist2  = _mm256_add_ps(_dx2dy2, _dz2);

                        // Mask particles within smoothing radius
                        _mask = _mm256_cmp_ps(_dist2, _smoothing_radius2, _CMP_LT_OS);

                        // Compute density for particles within the smoothing radius
                        //pi.density += pj.mass * settings.poly6 * std::pow(settings.smoothing_radius2 - distance2, 3.0f);
                        _radius_diff       = _mm256_sub_ps(_smoothing_radius2, _dist2);
                        _radius_diff2      = _mm256_mul_ps(_radius_diff, _radius_diff);
                        _radius_diff3      = _mm256_mul_ps(_radius_diff2, _radius_diff);
                        _poly6_radius_diff = _mm256_mul_ps(_poly6, _radius_diff3);
                        _density           = _mm256_mul_ps(_mass, _poly6_radius_diff);
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
        pi.pressure = settings.gas_constant * (pi.density - settings.rest_density);
    }
}

void SPH::ComputePressureViscousSurfaceTensionForce()
{
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        pi.acceleration = {};
        vf3 pressure_force = {};
        vf3 viscous_force = {};
        vf3 surface_tension_force = {};
        vf3 normal = {};
        vf3 gravity_force = pi.mass * settings.gravity;

        vi3 cell_origin = cell(pi, settings.smoothing_radius);

        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++)
            {
                for (s32 z = -1; z <= 1; z++)
                {
                    u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                    const auto& neighbors = particle_table[cell_hash];
                    // TODO: Intrinsics
                    for (s32 j = 0; j < neighbors.size(); j++)
                    {
                        Particle& pj = fluid_particles[neighbors[j]];

                        if (&pi == &pj || pj.hash != cell_hash)
                            continue;

                        vf3 difference = pj.position - pi.position;
                        f32 distance2 = glm::length2(difference);
                        if (distance2 < settings.smoothing_radius2)
                        {
                            f32 distance = std::sqrt(distance2);
                            vf3 direction = glm::normalize(difference);

                            pressure_force += -direction * pj.mass * (pi.pressure + pj.pressure) / (2.0f * pj.density) * settings.spiky_grad * std::pow(settings.smoothing_radius - distance, 3.0f);
                            viscous_force += pj.mass * settings.viscosity * (pj.velocity - pi.velocity) / pj.density * settings.spiky_laplacian * (settings.smoothing_radius - distance);
                            normal += direction * pj.mass / pj.density * settings.spiky_grad * std::pow(settings.smoothing_radius - distance, 3.0f);
                        }
                    }
                }
            }
        }

        f32 curvature = -glm::length(normal);
        if (curvature > 0.0f)
            surface_tension_force = -settings.surface_tension_constant * curvature * settings.spiky_laplacian * glm::normalize(normal);

        // REVISE: Do we add the surface tension force before or after the density division?
        pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;
    }
}

void SPH::ComputePressureViscousSurfaceTensionForceIntrinsics()
{
    __m256 _gravity_x = _mm256_set1_ps(settings.gravity.x);
    __m256 _gravity_y = _mm256_set1_ps(settings.gravity.y);
    __m256 _gravity_z = _mm256_set1_ps(settings.gravity.z);

    __m256 _smoothing_radius = _mm256_set1_ps(settings.smoothing_radius);
    __m256 _smoothing_radius2 = _mm256_set1_ps(settings.smoothing_radius2);
    __m256 _spiky_grad = _mm256_set1_ps(settings.spiky_grad);
    __m256 _spiky_laplacian = _mm256_set1_ps(settings.spiky_laplacian);
    __m256 _viscosity = _mm256_set1_ps(settings.viscosity);
    __m256 _two = _mm256_set1_ps(2.0f);

    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        __m256 _acc_x = _mm256_setzero_ps();
        __m256 _acc_y = _mm256_setzero_ps();
        __m256 _acc_z = _mm256_setzero_ps();

        __m256 _press_x = _mm256_setzero_ps();
        __m256 _press_y = _mm256_setzero_ps();
        __m256 _press_z = _mm256_setzero_ps();

        __m256 _visc_x = _mm256_setzero_ps();
        __m256 _visc_y = _mm256_setzero_ps();
        __m256 _visc_z = _mm256_setzero_ps();

        __m256 _surf_x = _mm256_setzero_ps();
        __m256 _surf_y = _mm256_setzero_ps();
        __m256 _surf_z = _mm256_setzero_ps();

        __m256 _norm_x = _mm256_setzero_ps();
        __m256 _norm_y = _mm256_setzero_ps();
        __m256 _norm_z = _mm256_setzero_ps();

        __m256 _gravity_force_x = _mm256_mul_ps(_mm256_set1_ps(pi.mass), _gravity_x);
        __m256 _gravity_force_y = _mm256_mul_ps(_mm256_set1_ps(pi.mass), _gravity_y);
        __m256 _gravity_force_z = _mm256_mul_ps(_mm256_set1_ps(pi.mass), _gravity_z);

        vi3 cell_origin = cell(pi, settings.smoothing_radius);

        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++)
            {
                for (s32 z = -1; z <= 1; z++)
                {
                    u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                    const auto& neighbors = particle_table[cell_hash];

                    for (s32 j = 0; j < neighbors.size(); j += 8)
                    {
                        // Load neighboring particle data
                        __m256 _px = _mm256_loadu_ps(&fluid_particles[neighbors[j]].position.x);
                        __m256 _py = _mm256_loadu_ps(&fluid_particles[neighbors[j]].position.y);
                        __m256 _pz = _mm256_loadu_ps(&fluid_particles[neighbors[j]].position.z);

                        __m256 _vx = _mm256_loadu_ps(&fluid_particles[neighbors[j]].velocity.x);
                        __m256 _vy = _mm256_loadu_ps(&fluid_particles[neighbors[j]].velocity.y);
                        __m256 _vz = _mm256_loadu_ps(&fluid_particles[neighbors[j]].velocity.z);

                        __m256 _pj_mass = _mm256_loadu_ps(&fluid_particles[neighbors[j]].mass);
                        __m256 _pj_pressure = _mm256_loadu_ps(&fluid_particles[neighbors[j]].pressure);
                        __m256 _pj_density = _mm256_loadu_ps(&fluid_particles[neighbors[j]].density);

                        // Compute difference
                        __m256 _dx = _mm256_sub_ps(_px, _mm256_set1_ps(pi.position.x));
                        __m256 _dy = _mm256_sub_ps(_py, _mm256_set1_ps(pi.position.y));
                        __m256 _dz = _mm256_sub_ps(_pz, _mm256_set1_ps(pi.position.z));

                        // Compute distance squared
                        __m256 _distance2 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_dx, _dx), _mm256_mul_ps(_dy, _dy)), _mm256_mul_ps(_dz, _dz));

                        // Mask for particles within smoothing radius
                        __m256 _mask = _mm256_cmp_ps(_distance2, _smoothing_radius2, _CMP_LT_OS);

                        // Distance (with sqrt)
                        __m256 _distance = _mm256_sqrt_ps(_distance2);

                        // Direction
                        __m256 inv_distance = _mm256_div_ps(_mm256_set1_ps(1.0f), _distance);
                        __m256 dir_x = _mm256_mul_ps(_dx, inv_distance);
                        __m256 dir_y = _mm256_mul_ps(_dy, inv_distance);
                        __m256 dir_z = _mm256_mul_ps(_dz, inv_distance);

                        // Pressure force
                        __m256 pressure_factor = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_sub_ps(_smoothing_radius, _distance), _mm256_sub_ps(_smoothing_radius, _distance)), _mm256_sub_ps(_smoothing_radius, _distance)), _spiky_grad);
                        pressure_factor = _mm256_mul_ps(pressure_factor, _mm256_div_ps(_mm256_add_ps(_mm256_set1_ps(pi.pressure), _pj_pressure), _mm256_mul_ps(_two, _pj_density)));
                        pressure_factor = _mm256_mul_ps(pressure_factor, _pj_mass);

                        pressure_factor = _mm256_blendv_ps(_mm256_setzero_ps(), pressure_factor, _mask);

                        _press_x = _mm256_fmadd_ps(pressure_factor, dir_x, _press_x);
                        _press_y = _mm256_fmadd_ps(pressure_factor, dir_y, _press_y);
                        _press_z = _mm256_fmadd_ps(pressure_factor, dir_z, _press_z);

                        // Viscous force
                        __m256 visc_factor = _mm256_mul_ps(_mm256_mul_ps(_pj_mass, _viscosity), _mm256_mul_ps(_mm256_sub_ps(_smoothing_radius, _distance), _spiky_laplacian));
                        visc_factor = _mm256_div_ps(visc_factor, _pj_density);

                        visc_factor = _mm256_blendv_ps(_mm256_setzero_ps(), visc_factor, _mask);

                        _visc_x = _mm256_fmadd_ps(visc_factor, _mm256_sub_ps(_vx, _mm256_set1_ps(pi.velocity.x)), _visc_x);
                        _visc_y = _mm256_fmadd_ps(visc_factor, _mm256_sub_ps(_vy, _mm256_set1_ps(pi.velocity.y)), _visc_y);
                        _visc_z = _mm256_fmadd_ps(visc_factor, _mm256_sub_ps(_vz, _mm256_set1_ps(pi.velocity.z)), _visc_z);

                        // Normal (for surface tension)
                        _norm_x = _mm256_fmadd_ps(pressure_factor, dir_x, _norm_x);
                        _norm_y = _mm256_fmadd_ps(pressure_factor, dir_y, _norm_y);
                        _norm_z = _mm256_fmadd_ps(pressure_factor, dir_z, _norm_z);
                    }
                }
            }
        }

        // Final forces and acceleration
        __m256 curvature = _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_norm_x, _norm_x), _mm256_mul_ps(_norm_y, _norm_y)), _mm256_mul_ps(_norm_z, _norm_z)));
        curvature = _mm256_sub_ps(_mm256_set1_ps(0.0f), curvature);

        if (_mm256_movemask_ps(_mm256_cmp_ps(curvature, _mm256_set1_ps(0.0f), _CMP_GT_OS)))
        {
            _surf_x = _mm256_mul_ps(curvature, _norm_x);
            _surf_y = _mm256_mul_ps(curvature, _norm_y);
            _surf_z = _mm256_mul_ps(curvature, _norm_z);
        }

        // Update particle acceleration
        _acc_x = _mm256_add_ps(_acc_x, _mm256_div_ps(_mm256_add_ps(_mm256_add_ps(_press_x, _visc_x), _surf_x), _mm256_set1_ps(pi.density)));
        _acc_y = _mm256_add_ps(_acc_y, _mm256_div_ps(_mm256_add_ps(_mm256_add_ps(_press_y, _visc_y), _surf_y), _mm256_set1_ps(pi.density)));
        _acc_z = _mm256_add_ps(_acc_z, _mm256_div_ps(_mm256_add_ps(_mm256_add_ps(_press_z, _visc_z), _surf_z), _mm256_set1_ps(pi.density)));

        _acc_x = _mm256_add_ps(_acc_x, _gravity_force_x);
        _acc_y = _mm256_add_ps(_acc_y, _gravity_force_y);
        _acc_z = _mm256_add_ps(_acc_z, _gravity_force_z);

        // Store results back to particle
        _mm256_storeu_ps(&pi.acceleration.x, _acc_x);
        _mm256_storeu_ps(&pi.acceleration.y, _acc_y);
        _mm256_storeu_ps(&pi.acceleration.z, _acc_z);
    }

}

void SPH::Integrate(f32 dt)
{
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        // Update Velocity
        pi.velocity += pi.acceleration * dt;

        // Update Position
        pi.position += pi.velocity * dt;
    }
}

void SPH::ComputeBoundaryCondition()
{
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        // Particle Collision
        vi3 cell_origin = cell(pi, settings.smoothing_radius);

        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++)
            {
                for (s32 z = -1; z <= 1; z++)
                {
                    u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                    const auto& neighbors = particle_table[cell_hash];
                    //TODO: Intrinsics
                    for (s32 j = 0; j < neighbors.size(); j++)
                    {
                        Particle& pj = fluid_particles[neighbors[j]];

                        if (&pi == &pj || pj.hash != cell_hash)
                            continue;

                        vf3 difference = pj.position - pi.position;
                        f32 distance2 = glm::length2(difference);
                        f32 distance = std::sqrt(distance2);

                        if (distance < 0.0f)
                        {
                            if (distance < 1e-6f)
                                distance = 1e-6f;

                            vf3 direction = difference / distance;
                            f32 overlap = settings.smoothing_radius - distance;

                            // Move particles apart to resolve overlap
                            pi.position -= direction * overlap * 0.5f;
                            pj.position += direction * overlap * 0.5f;

                            // Compute the relative velocity in the direction of the collision
                            vf3 relative_velocity = pi.velocity - pj.velocity;
                            f32 velocity_along_normal = glm::dot(relative_velocity, direction);

                            // Calculate the impulse to be applied
                            f32 combined_mass = pi.mass + pj.mass;
                            f32 impulse = (2.0f * velocity_along_normal) / combined_mass;

                            // Apply the impulse to the velocities of both particles
                            pi.velocity -= impulse * (pj.mass / combined_mass) * direction;
                            pj.velocity += impulse * (pi.mass / combined_mass) * direction;
                        }
                    }
                }
            }
        }

        // Boundary Collision
        if (pi.position.x - settings.boundary_epsilon < settings.boundary_min.x)
        {
            pi.velocity.x *= settings.boundary_damping;
            pi.position.x = settings.boundary_min.x + settings.boundary_epsilon;
        }
        if (pi.position.x + settings.boundary_epsilon > settings.boundary_max.x)
        {
            pi.velocity.x *= settings.boundary_damping;
            pi.position.x = settings.boundary_max.x - settings.boundary_epsilon;
        }

        if (pi.position.y - settings.boundary_epsilon < settings.boundary_min.y)
        {
            pi.velocity.y *= settings.boundary_damping;
            pi.position.y = settings.boundary_min.y + settings.boundary_epsilon;
        }
        if (pi.position.y + settings.boundary_epsilon > settings.boundary_max.y)
        {
            pi.velocity.y *= settings.boundary_damping;
            pi.position.y = settings.boundary_max.y - settings.boundary_epsilon;
        }

        if (pi.position.z - settings.boundary_epsilon < settings.boundary_min.z)
        {
            pi.velocity.z *= settings.boundary_damping;
            pi.position.z = settings.boundary_min.z + settings.boundary_epsilon;
        }
        if (pi.position.z + settings.boundary_epsilon > settings.boundary_max.z)
        {
            pi.velocity.z *= settings.boundary_damping;
            pi.position.z = settings.boundary_max.z - settings.boundary_epsilon;
        }
    }
}

std::vector<Particle>& SPH::GetParticles() { return fluid_particles; }

void SPH::CreateNeighborTable()
{
    // Create Neighbor Table
    particle_table.clear();
    particle_table.resize(TABLE_SIZE);
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        u32 current_hash = fluid_particles[i].hash;
        particle_table[current_hash].push_back(i);
    }
}

// Hash map
vi3 SPH::cell(const Particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
u32 SPH::hash(const vi3& cell) { return (static_cast<u32>(cell.x * 92837111) ^ static_cast<u32>(cell.y * 689287499) ^ static_cast<u32>(cell.z * 283823481)) % TABLE_SIZE; }
