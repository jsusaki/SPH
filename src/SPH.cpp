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
    //ComputePressureViscousSurfaceTensionForce();

    ComputeDensityPressureIntrinsics();
    ComputePressureViscousSurfaceTensionForceIntrinsics();

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
    f32 pj_px[8] = { 0.0f }, pj_py[8] = { 0.0f }, pj_pz[8] = { 0.0f };

    // 64-bit float registers
    // Positions
    __m256 _mass, _density;
    __m256 _pi_px, _pi_py, _pi_pz;
    __m256 _pj_px, _pj_py, _pj_pz;

    // Squared distance
    __m256 _dx, _dy, _dz, _dx2, _dy2, _dz2, _dx2dy2;
    __m256 _dist2;
    // if condition with bitwise AND operation
    __m256 _mask;
    // Kernels
    __m256 _poly6, _smoothing_radius2;
    __m256 _radius_diff, _radius_diff2, _radius_diff3, _poly6_radius_diff;

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
        _pi_px = _mm256_set1_ps(pi.position.x);
        _pi_py = _mm256_set1_ps(pi.position.y);
        _pi_pz = _mm256_set1_ps(pi.position.z);

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
                        for (s32 k = 0; k < 8; k++)
                        {
                            if (j + k < neighbors.size())
                            {
                                const Particle& pj = fluid_particles[neighbors[j + k]];
                                pj_px[k] = pj.position.x;
                                pj_py[k] = pj.position.y;
                                pj_pz[k] = pj.position.z;
                            }
                        }

                        // Load neighbor positions
                        _pj_px = _mm256_loadu_ps(pj_px);
                        _pj_py = _mm256_loadu_ps(pj_py);
                        _pj_pz = _mm256_loadu_ps(pj_pz);

                        // Compute squared distance
                        _dx     = _mm256_sub_ps(_pi_px, _pj_px);
                        _dy     = _mm256_sub_ps(_pi_py, _pj_py);
                        _dz     = _mm256_sub_ps(_pi_pz, _pj_pz);
                        _dx2    = _mm256_mul_ps(_dx, _dx);
                        _dy2    = _mm256_mul_ps(_dy, _dy);
                        _dz2    = _mm256_mul_ps(_dz, _dz);
                        _dx2dy2 = _mm256_add_ps(_dx2, _dy2);
                        _dist2  = _mm256_add_ps(_dx2dy2, _dz2);

                        // Mask particles within smoothing radius
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
                        f32 density[8] = { 0.0f };
                        _mm256_storeu_ps(density, _density);

                        f32 density_sum = 0.0f;
                        for (s32 j = 0; j < 8; j++)
                        {
                            if (!std::isnan(density[j]))
                                density_sum += density[j];
                        }
                        pi.density += density_sum;
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

                            pressure_force += -direction * pj.mass * ((pi.pressure + pj.pressure) / (2.0f * pj.density)) * settings.spiky_grad * std::pow(settings.smoothing_radius - distance, 3.0f);
                            viscous_force += pj.mass * settings.viscosity * ((pj.velocity - pi.velocity) / pj.density) * settings.spiky_laplacian * (settings.smoothing_radius - distance);
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
    f32 pj_px[8] = { 0.0f }, pj_py[8] = { 0.0f }, pj_pz[8] = { 0.0f };
    f32 pj_vx[8] = { 0.0f }, pj_vy[8] = { 0.0f }, pj_vz[8] = { 0.0f };
    f32 pj_density[8] = { 0.0f }, pj_pressure[8] = { 0.0f };

    __m256 _viscosity, _mass;
    __m256 _gravity_x, _gravity_y, _gravity_z;

    // pi data
    __m256 _pi_px, _pi_py, _pi_pz;
    __m256 _pi_vx, _pi_vy, _pi_vz;
    __m256 _pi_pressure, _pi_density;

    // pj data
    __m256 _pj_px, _pj_py, _pj_pz;
    __m256 _pj_vx, _pj_vy, _pj_vz;
    __m256 _pj_density, _pj_pressure;

    // Squared distance
    __m256 _dx, _dy, _dz, _dx2, _dy2, _dz2, _dx2dy2;
    __m256 _dist2, _dist;
    __m256 _inv_dist, _dir_x, _dir_y, _dir_z;

    // Kernels
    __m256 _smoothing_radius, _smoothing_radius2, _spiky_grad, _spiky_laplacian;
    __m256 _radius_diff, _radius_diff2, _radius_diff3;

    // Forces
    __m256 _acceleration_x, _acceleration_y, _acceleration_z;
    __m256 _gravity_force_x, _gravity_force_y, _gravity_force_z;

    // Pressure
    __m256 _pressure_component, _pressure_force_x, _pressure_force_y, _pressure_force_z;
    __m256 _pressure_density, _two_pj_density, _mass_pressure_density_spiky_grad, _pi_pj_pressure_add, _mass_pressure_density;
    
    // Viscosity
    __m256 _viscous_force, _viscous_force_x, _viscous_force_y, _viscous_force_z;
    __m256 _pj_mass_visc, _spiky_laplacian_pj_density, _spiky_laplacian_pj_density_radius_diff;
    __m256 _dvx, _dvy, _dvz, _dvx_density, _dvy_density, _dvz_density;

    // Surface Tension
    __m256 _curvature, _surface_tension_x, _surface_tension_y, _surface_tension_z;
    __m256 _normal_x, _normal_y, _normal_z;
    __m256 _norm_mag, _inv_norm_mag, _norm_x_unit, _norm_y_unit, _norm_z_unit;

    // Constants
    __m256 _zero, _one, _two, _mask1, _mask2;

    __m256 _min, _max;

    _smoothing_radius  = _mm256_set1_ps(settings.smoothing_radius);
    _smoothing_radius2 = _mm256_set1_ps(settings.smoothing_radius2);
    _spiky_grad        = _mm256_set1_ps(settings.spiky_grad);
    _spiky_laplacian   = _mm256_set1_ps(settings.spiky_laplacian);

    _gravity_x = _mm256_set1_ps(settings.gravity.x);
    _gravity_y = _mm256_set1_ps(settings.gravity.y);
    _gravity_z = _mm256_set1_ps(settings.gravity.z);
    _viscosity = _mm256_set1_ps(settings.viscosity);
    _mass      = _mm256_set1_ps(settings.mass);

    _zero = _mm256_set1_ps(0.0f);
    _one  = _mm256_set1_ps(1.0f);
    _two  = _mm256_set1_ps(2.0f);

    _min = _mm256_set1_ps(std::numeric_limits<f32>::min());
    _max = _mm256_set1_ps(std::numeric_limits<f32>::max());

    //_min = _mm256_set1_ps(-10.0f);
    //_max = _mm256_set1_ps(10.0f);
    auto clamp = [](__m256 a, __m256 min, __m256 max) { return _mm256_max_ps(min, _mm256_min_ps(a, max)); };

    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        // Set own positions
        _pi_px = _mm256_set1_ps(pi.position.x);
        _pi_py = _mm256_set1_ps(pi.position.y);
        _pi_pz = _mm256_set1_ps(pi.position.z);

        _pi_vx = _mm256_set1_ps(pi.velocity.x);
        _pi_vy = _mm256_set1_ps(pi.velocity.y);
        _pi_vz = _mm256_set1_ps(pi.velocity.z);

        _pi_density  = _mm256_set1_ps(pi.density);
        _pi_pressure = _mm256_set1_ps(pi.pressure);

        // Reset acceleration
        _acceleration_x = _mm256_setzero_ps();
        _acceleration_y = _mm256_setzero_ps();
        _acceleration_z = _mm256_setzero_ps();
        
        // Reset pressure
        _pressure_component   = _mm256_setzero_ps();
        _pressure_force_x = _mm256_setzero_ps();
        _pressure_force_y = _mm256_setzero_ps();
        _pressure_force_z = _mm256_setzero_ps();

        // Reset viscous
        _viscous_force   = _mm256_setzero_ps();
        _viscous_force_x = _mm256_setzero_ps();
        _viscous_force_y = _mm256_setzero_ps();
        _viscous_force_z = _mm256_setzero_ps();

        // Reset surface tension
        _curvature         = _mm256_setzero_ps();
        _surface_tension_x = _mm256_setzero_ps();
        _surface_tension_y = _mm256_setzero_ps();
        _surface_tension_z = _mm256_setzero_ps();

        _normal_x = _mm256_setzero_ps();
        _normal_y = _mm256_setzero_ps();
        _normal_z = _mm256_setzero_ps();
        
        // Compute gravity
        _gravity_force_x = _mm256_mul_ps(_mass, _gravity_x);
        _gravity_force_y = _mm256_mul_ps(_mass, _gravity_y);
        _gravity_force_z = _mm256_mul_ps(_mass, _gravity_z);

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
                        for (s32 k = 0; k < 8; k++)
                        {
                            if (j + k < neighbors.size())
                            {
                                const Particle& pj = fluid_particles[neighbors[j + k]];
                                pj_px[k] = pj.position.x;
                                pj_py[k] = pj.position.y;
                                pj_pz[k] = pj.position.z;
                                pj_vx[k] = pj.velocity.x;
                                pj_vy[k] = pj.velocity.y;
                                pj_vz[k] = pj.velocity.z;
                                pj_density[k]  = pj.density;
                                pj_pressure[k] = pj.pressure;
                            }
                        }

                        // Load neighboring particle data
                        _pj_px = _mm256_loadu_ps(pj_px);
                        _pj_py = _mm256_loadu_ps(pj_py);
                        _pj_pz = _mm256_loadu_ps(pj_pz);

                        _pj_vx = _mm256_loadu_ps(pj_vx);
                        _pj_vy = _mm256_loadu_ps(pj_vy);
                        _pj_vz = _mm256_loadu_ps(pj_vz);

                        _pj_density  = _mm256_loadu_ps(pj_density);
                        _pj_pressure = _mm256_loadu_ps(pj_pressure);

                        // Compute squared distance
                        _dx     = _mm256_sub_ps(_pj_px, _pi_px);
                        _dy     = _mm256_sub_ps(_pj_py, _pi_py);
                        _dz     = _mm256_sub_ps(_pj_pz, _pi_pz);
                        _dx2    = _mm256_mul_ps(_dx, _dx);
                        _dy2    = _mm256_mul_ps(_dy, _dy);
                        _dz2    = _mm256_mul_ps(_dz, _dz);
                        _dx2dy2 = _mm256_add_ps(_dx2, _dy2);
                        _dist2  = _mm256_add_ps(_dx2dy2, _dz2);

                        // Mask particles within smoothing radius
                        _mask1 = _mm256_cmp_ps(_dist2, _smoothing_radius2, _CMP_LT_OS);
                        
                        // Distance
                        //f32 distance = sqrt(distance2);
                        _dist = _mm256_sqrt_ps(_dist2);
                        _inv_dist = _mm256_div_ps(_one, _dist);

                        // Pressure force
                        // pressure_force += -direction * pj.mass * ((pi.pressure + pj.pressure) / (2.0f * pj.density)) * spiky_grad * pow(smoothing_radius - distance, 3.0f);
                       
                        // pow(smoothing_radius - distance, 3.0f)
                        _radius_diff  = _mm256_sub_ps(_smoothing_radius, _dist);
                        _radius_diff2 = _mm256_mul_ps(_radius_diff, _radius_diff);
                        _radius_diff3 = _mm256_mul_ps(_radius_diff2, _radius_diff);

                        // (pi.pressure + pj.pressure)
                        _pi_pj_pressure_add = _mm256_add_ps(_pi_pressure, _pj_pressure);
                        // (2.0f * pj.density)
                        _two_pj_density = _mm256_mul_ps(_two, _pj_density);
                        // (pi.pressure + pj.pressure) / (2.0f * pj.density) 
                        _pressure_density = _mm256_div_ps(_pi_pj_pressure_add, _two_pj_density);
                        // pj.mass * ((pi.pressure + pj.pressure)  / (2.0f * pj.density)) 
                        _mass_pressure_density = _mm256_mul_ps(_mass, _pressure_density);
                        // pj.mass * ((pi.pressure + pj.pressure)  / (2.0f * pj.density)) * spiky_grad
                        _mass_pressure_density_spiky_grad = _mm256_mul_ps(_mass_pressure_density, _spiky_grad);
                        // pj.mass * ((pi.pressure + pj.pressure)  / (2.0f * pj.density)) * spiky_grad * pow(smoothing_radius - distance, 3.0f)
                        _pressure_component = _mm256_mul_ps(_mass_pressure_density_spiky_grad, _radius_diff3);
                        // if dist2 < radius2
                        _pressure_component = _mm256_blendv_ps(_zero, _pressure_component, _mask1);

                        //vf3 direction = glm::normalize(difference);
                        _dir_x = _mm256_mul_ps(_dx, _inv_dist);
                        _dir_y = _mm256_mul_ps(_dy, _inv_dist);
                        _dir_z = _mm256_mul_ps(_dz, _inv_dist);

                        // pressure_force += -direction * pj.mass * ((pi.pressure + pj.pressure) / (2.0f * pj.density)) * spiky_grad * pow(smoothing_radius - distance, 3.0f);
                        _pressure_force_x = _mm256_fmadd_ps(_pressure_component, _dir_x, _pressure_force_x);
                        _pressure_force_y = _mm256_fmadd_ps(_pressure_component, _dir_y, _pressure_force_y);
                        _pressure_force_z = _mm256_fmadd_ps(_pressure_component, _dir_z, _pressure_force_z);

                        // Viscous force
                        // viscous_force += pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density) * spiky_laplacian * (smoothing_radius - distance);

                        // (pj.velocity - pi.velocity)
                        _dvx = _mm256_sub_ps(_pj_vx, _pi_vx);
                        _dvy = _mm256_sub_ps(_pj_vy, _pi_vy);
                        _dvz = _mm256_sub_ps(_pj_vz, _pi_vz);

                        // ((pj.velocity - pi.velocity) / pj.density)
                        _dvx_density = _mm256_div_ps(_dvx, _pj_density);
                        _dvy_density = _mm256_div_ps(_dvy, _pj_density);
                        _dvz_density = _mm256_div_ps(_dvz, _pj_density);

                        // viscosity * ((pj.velocity - pi.velocity) / pj.density)
                        __m256 _viscous_comp_x = _mm256_mul_ps(_dvx_density, _viscosity);
                        __m256 _viscous_comp_y = _mm256_mul_ps(_dvy_density, _viscosity);
                        __m256 _viscous_comp_z = _mm256_mul_ps(_dvz_density, _viscosity);

                        _viscous_comp_x = _mm256_mul_ps(_viscous_comp_x, _mass);
                        _viscous_comp_y = _mm256_mul_ps(_viscous_comp_y, _mass);
                        _viscous_comp_z = _mm256_mul_ps(_viscous_comp_z, _mass);

                        _viscous_comp_x = _mm256_mul_ps(_viscous_comp_x, _spiky_laplacian);
                        _viscous_comp_y = _mm256_mul_ps(_viscous_comp_y, _spiky_laplacian);
                        _viscous_comp_z = _mm256_mul_ps(_viscous_comp_z, _spiky_laplacian);

                        _viscous_comp_x = _mm256_mul_ps(_viscous_comp_x, _radius_diff);
                        _viscous_comp_y = _mm256_mul_ps(_viscous_comp_y, _radius_diff);
                        _viscous_comp_z = _mm256_mul_ps(_viscous_comp_z, _radius_diff);

                        // Apply the mask to the viscous force components
                        _viscous_comp_x = _mm256_blendv_ps(_zero, _viscous_comp_x, _mask1);
                        _viscous_comp_y = _mm256_blendv_ps(_zero, _viscous_comp_y, _mask1);
                        _viscous_comp_z = _mm256_blendv_ps(_zero, _viscous_comp_z, _mask1);

                        // viscous_force += pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density) * spiky_laplacian * (smoothing_radius - distance);
                        _viscous_force_x = _mm256_add_ps(_viscous_force_x, _viscous_comp_x);
                        _viscous_force_y = _mm256_add_ps(_viscous_force_y, _viscous_comp_y);
                        _viscous_force_z = _mm256_add_ps(_viscous_force_z, _viscous_comp_z);

                        // Normal
                        _normal_x = _mm256_fmadd_ps(_pressure_component, _dir_x, _normal_x);
                        _normal_y = _mm256_fmadd_ps(_pressure_component, _dir_y, _normal_y);
                        _normal_z = _mm256_fmadd_ps(_pressure_component, _dir_z, _normal_z);
                    }
                }
            }
        }

        //f32 curvature = -glm::length(normal);
        //if (curvature > 0.0f)
            //surface_tension_force = -settings.surface_tension_constant * curvature * settings.spiky_laplacian * glm::normalize(normal);
        
        // Recalculate norm vector using density gradient for surface tension

        _norm_mag = _mm256_sqrt_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(_normal_x, _normal_x), 
                    _mm256_mul_ps(_normal_y, _normal_y)), 
                _mm256_mul_ps(_normal_z, _normal_z)
            ));
        
        __m256 _epsilon = _mm256_set1_ps(1e-6f);
        _norm_mag = _mm256_max_ps(_norm_mag, _epsilon);
        _inv_norm_mag = _mm256_div_ps(_one, _norm_mag);

        // Normalize the normal vector
        _inv_norm_mag = _mm256_div_ps(_one, _norm_mag);
        _norm_x_unit = _mm256_mul_ps(_normal_x, _inv_norm_mag);
        _norm_y_unit = _mm256_mul_ps(_normal_y, _inv_norm_mag);
        _norm_z_unit = _mm256_mul_ps(_normal_z, _inv_norm_mag);

        // Surface tension curvature
        _curvature  = _mm256_sub_ps(_zero, _norm_mag); // Negative curvature indicates a concave surface
        _mask2 = _mm256_cmp_ps(_curvature, _zero, _CMP_GT_OS);

        // Surface tension forces
        _surface_tension_x = _mm256_blendv_ps(_surface_tension_x, _mm256_mul_ps(_curvature, _norm_x_unit), _mask2);
        _surface_tension_y = _mm256_blendv_ps(_surface_tension_y, _mm256_mul_ps(_curvature, _norm_y_unit), _mask2);
        _surface_tension_z = _mm256_blendv_ps(_surface_tension_z, _mm256_mul_ps(_curvature, _norm_z_unit), _mask2);

        // Update particle acceleration
        // pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;

        // Update Pressure Force
        // acceleration += pressure_force
        _acceleration_x = _mm256_add_ps(_acceleration_x, _pressure_force_x);
        _acceleration_y = _mm256_add_ps(_acceleration_y, _pressure_force_y);
        _acceleration_z = _mm256_add_ps(_acceleration_z, _pressure_force_z);

        // acceleration += viscous_force
        _acceleration_x = _mm256_add_ps(_acceleration_x, _viscous_force_x);
        _acceleration_y = _mm256_add_ps(_acceleration_y, _viscous_force_y);
        _acceleration_z = _mm256_add_ps(_acceleration_z, _viscous_force_z);

        // acceleration += surface_tension
        _acceleration_x = _mm256_add_ps(_acceleration_x, _surface_tension_x);
        _acceleration_y = _mm256_add_ps(_acceleration_y, _surface_tension_y);
        _acceleration_z = _mm256_add_ps(_acceleration_z, _surface_tension_z);

        // acceleration /= pi.density
        _acceleration_x = _mm256_div_ps(_acceleration_x, _pi_density);
        _acceleration_y = _mm256_div_ps(_acceleration_y, _pi_density);
        _acceleration_z = _mm256_div_ps(_acceleration_z, _pi_density);

        // acceleration += gravity_force
        _acceleration_x = _mm256_add_ps(_acceleration_x, _gravity_force_x);
        _acceleration_y = _mm256_add_ps(_acceleration_y, _gravity_force_y);
        _acceleration_z = _mm256_add_ps(_acceleration_z, _gravity_force_z);

        // Store the results from SIMD registers to arrays
        f32 acceleration_x[8] = { 0.0f };
        f32 acceleration_y[8] = { 0.0f };
        f32 acceleration_z[8] = { 0.0f };
        _mm256_storeu_ps(acceleration_x, _acceleration_x);
        _mm256_storeu_ps(acceleration_y, _acceleration_y);
        _mm256_storeu_ps(acceleration_z, _acceleration_z);

        // Loop through the array to store the results back to pi.acceleration for each particle
        f32 acc_x_sum = 0.0f, acc_y_sum = 0.0f, acc_z_sum = 0.0f;
        for (s32 j = 0; j < 8; j++)
        {
            if (!std::isnan(acceleration_x[j])) acc_x_sum += acceleration_x[j];
            if (!std::isnan(acceleration_y[j])) acc_y_sum += acceleration_y[j];
            if (!std::isnan(acceleration_z[j])) acc_z_sum += acceleration_z[j];
        }

        pi.acceleration.x += acc_x_sum;
        pi.acceleration.y += acc_y_sum;
        pi.acceleration.z += acc_z_sum;
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
    u32 prev_hash = NO_PARTICLE;
        for (s32 i = 0; i < fluid_particles.size(); i++)
        {
            u32 current_hash = fluid_particles[i].hash;
            if (current_hash != prev_hash)
            {
                particle_table[current_hash].push_back(i);
                prev_hash = current_hash;
            }
        }
}

// Hash map
vi3 SPH::cell(const Particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
u32 SPH::hash(const vi3& cell) { return (static_cast<u32>(cell.x * 92837111) ^ static_cast<u32>(cell.y * 689287499) ^ static_cast<u32>(cell.z * 283823481)) % TABLE_SIZE; }
