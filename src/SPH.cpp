#include "SPH.h"

// Correctness
    // TODO: Review Particle Collision Physics
    // TODO: Debugging Info (support velocity field color, support radius color, collision contact point color)
        // DONE: Velocity Field Color
    // TODO: Density, Pressure Intrinsics
    // TODO: Pressure Force + Viscosity Force + Surface Tension Intrinsics

// Graphics
    // DONE: Velocity-based Particle Color
    // TODO: Ability to change color palette
    // TODO: Density-based Particle Radius Color
    // TODO: Marching Cube Algorithm for surface polygon mesh generation
    // TODO: Fluid Simulation Scenes: Dam Break, etc.
    // TODO: Teapot Scene, Bruce Lee, Be formless, shapeless, like water be water my friend
    // TODO: Light Shader

// Optimization
    // TODO: Boundary Condition: Collision SIMD
    // TODO: Spatial Hash Grid Search
    // TODO: SPH in Cuda

// Refactor
    // TODO: Particle Collision code transform into optimal form
    // TODO: Shader Management

// Features
    // TODO: Real-time Particle Repositioning
    // TODO: Mouse Right Button Force interaction
    // TODO: Interpolate Camera Positions
    // TODO: Physics-Informed Neural Network in C++ and CUDA

SPH::SPH() {}

SPH::SPH(const SPHSettings& s) { Init(s); }

void SPH::Init(const SPHSettings& s)
{
    settings = s;

    // TODO: real-time particle repositioning

    // Initialize Fluid Particles
    fluid_particles.clear();
    fluid_particles.reserve(settings.n_particles);

    // Particle spacing
    // Cube root to define the number of particles
    s32 cube_size      = static_cast<s32>(std::cbrt(settings.n_particles));
    f32 offset         = settings.radius + s.particle_spacing;
    f32 half_cube_size = (cube_size - 1) * offset / 2.0f;

    for (s32 x = 0; x < cube_size; x++)
    {
        for (s32 y = 0; y < cube_size; y++)
        {
            for (s32 z = 0; z < cube_size; z++)
            {
                if (fluid_particles.size() >= settings.n_particles)
                    break;

                Particle p;
                p.position = {
                    x * offset - half_cube_size + s.cube_offset.x,
                    y * offset - half_cube_size + s.cube_offset.y,
                    z * offset - half_cube_size + s.cube_offset.z,
                };
                p.velocity     = { 0.0f, 0.0f, 0.0f };
                p.acceleration = { 0.0f, 0.0f, 0.0f };
                p.mass         = s.mass;
                p.pressure     = 0.0f;
                p.density      = s.initial_particle_density;
                p.radius       = s.radius;
                p.support_radius = s.support_radius;
                p.viscosity    = s.viscosity;
                p.hash         = hash(cell(p, settings.support_radius));

                fluid_particles.push_back(p);
            }
        }
    }

    // Debug
    std::printf("grid_size=%d\n",   cube_size);
    std::printf("n_particles=%d\n", fluid_particles.size());

    // Sort Particles
    std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

    CreateNeighborTable();

    // TODO: Initialize Boundary Particles (obstacle)

    // Reset Statistics
    ComputeStatistics();
}


void SPH::Simulate(const SPHSettings& s)
{
    settings = s;
    // Compute hash
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& p = fluid_particles[i];
        p.hash = hash(cell(p, settings.support_radius));
    }

    // Sort Particles
    std::sort(fluid_particles.begin(), fluid_particles.end(), [&](const Particle& i, const Particle& j) { return i.hash < j.hash; });

    // Create Neighbor Table
    CreateNeighborTable();

    //TODO: Precompute neighbor with search?

    // Enable AVX2
    if (avx2)
    {
        ComputeDensityPressureIntrinsics();
        ComputePressureViscousSurfaceTensionForceIntrinsics();
    }
    else
    {
        ComputeDensityPressure();
        ComputePressureViscousSurfaceTensionForce();
    }

    Integrate(settings.dt);

    ComputeBoundaryCondition();

    ComputeStatistics();
}

void SPH::ComputeDensityPressure()
{
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        pi.density  = settings.initial_particle_density;
        pi.pressure = 0.0f;

        // Spatial Hash Grid Search
        vi3 cell_origin = cell(pi, settings.support_radius);
        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++)
            {
                for (s32 z = -1; z <= 1; z++)
                {
                    u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                    const auto& neighbors = particle_table[cell_hash];

                    // For each neighbor, compute the density
                    for (s32 j = 0; j < neighbors.size(); j++)
                    {
                        Particle& pj = fluid_particles[neighbors[j]];

                        if (&pi == &pj || cell_hash != pj.hash)
                            continue;

                        // Compute density for particles within the support radius
                        f32 dist2 = glm::distance2(pj.position, pi.position);
                        if (dist2 <= settings.support_radius2)
                        {
                            pi.density += pj.mass * settings.poly6 * std::pow(settings.support_radius2 - dist2, 3.0f);
                        }
                    }
                }
            }
        }

        // Compute pressure based on density
        pi.pressure = settings.stiffness * (pi.density - settings.rest_density);
    }
}

void SPH::ComputeDensityPressureIntrinsics()
{
    // Particle data arrays
    f32 pj_px[8] = { 0.0f }, pj_py[8] = { 0.0f }, pj_pz[8] = { 0.0f };

    // 32-bit float registers
    __m256 _mass, _density;

    __m256 _pi_px, _pi_py, _pi_pz;
    __m256 _pj_px, _pj_py, _pj_pz;

    __m256 _dx, _dy, _dz, _dx2, _dy2, _dz2, _dx2dy2;
    __m256 _dist2;

    __m256 _mask;

    __m256 _poly6, _support_radius2;
    __m256 _dr, _dr2, _dr3, _da;

    // Set constants into AVX2 registers
    _mass            = _mm256_set1_ps(settings.mass);
    _poly6           = _mm256_set1_ps(settings.poly6);
    _support_radius2 = _mm256_set1_ps(settings.support_radius2);

    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        pi.density  = settings.rest_density;
        pi.pressure = 0.0f;

        // Set own positions
        _pi_px = _mm256_set1_ps(pi.position.x);
        _pi_py = _mm256_set1_ps(pi.position.y);
        _pi_pz = _mm256_set1_ps(pi.position.z);

        // Spatial Hash Grid Search
        vi3 cell_origin = cell(pi, settings.support_radius);
        for (s32 x = -1; x <= 1; x++)
        {
            for (s32 y = -1; y <= 1; y++)
            {
                for (s32 z = -1; z <= 1; z++)
                {
                    u32 cell_hash = hash(cell_origin + vi3(x, y, z));
                    const auto& neighbors = particle_table[cell_hash];

                    // For each neighbor, compute the density
                    for (s32 j = 0; j < neighbors.size(); j += 8)
                    {
                        // Convert neighbor positions into array for SIMD
                        // TODO: there is a better way to achieve this? Instead of Array of Structre, make Structure of Arrays?
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

                        // Mask particles within support_radius
                        _mask   = _mm256_cmp_ps(_dist2, _support_radius2, _CMP_LT_OS);

                        // TODO: REINPLEMENT DENSITY, HAS SOME BUGS

                        // Compute density for particles within the support radius
                        //pi.density += pj.mass * settings.poly6 * pow(support_radius2 - dist2, 3.0f);
                        _dr      = _mm256_sub_ps(_support_radius2, _dist2);
                        _dr2     = _mm256_mul_ps(_dr, _dr);
                        _dr3     = _mm256_mul_ps(_dr2, _dr);
                        _da      = _mm256_mul_ps(_poly6, _dr3);
                        _density = _mm256_mul_ps(_mass, _da);
                        
                        // Apply mask to zero out the contribution for particles outside the support_radius
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

        // Compute pressure based on density
        pi.pressure = settings.stiffness * (pi.density - settings.rest_density);
    }
}

void SPH::ComputePressureViscousSurfaceTensionForce()
{
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        pi.acceleration    = {};
        vf3 pressure_force = {};
        vf3 viscous_force  = {};

        vf3 surface_tension_force = {};
        f32 surface_tension_lap   = 0.0f;
        vf3 surface_tension_grad  = {};
        vf3 gravity_force = pi.mass * settings.gravity;

        // Spatial Hash Grid Search
        vi3 cell_origin = cell(pi, settings.support_radius);
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

                        vf3 diff  = pj.position - pi.position;
                        f32 dist2 = glm::length2(diff);
                        if (dist2 <= settings.support_radius2)
                        {
                            f32 dist = std::sqrt(dist2);

                            pressure_force       += pj.mass * diff * ((pi.pressure + pj.pressure) / (2.0f * pj.density)) * settings.spiky_grad * std::pow(settings.support_radius - dist, 3.0f) / dist;
                            viscous_force        += pj.mass * settings.viscosity * ((pj.velocity - pi.velocity) / pj.density) * settings.viscosity_laplacian * (settings.support_radius - dist);
                            
                            surface_tension_grad += (pj.mass / pj.density) * diff * settings.poly6_grad * std::pow(settings.support_radius2 - dist2, 2.0f);
                            surface_tension_lap  += (pj.mass / pj.density) * settings.poly6_grad * (settings.support_radius2 - dist2) * (3.0f * settings.support_radius2 - 7.0f * dist2);
                        }
                    }
                }
            }
        }

        f32 grad_length = glm::length(surface_tension_grad);
        f32 epsilon = 1e-6f;
        if (grad_length >= epsilon)
            surface_tension_force = -settings.surface_tension_constant * surface_tension_grad / (grad_length + epsilon) * surface_tension_lap;

        pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;
    }
}

void SPH::ComputePressureViscousSurfaceTensionForceIntrinsics()
{
    // Particle data arrays
    f32 pj_density[8] = { 0.0f }, pj_pressure[8] = { 0.0f };
    f32 pj_px[8] = { 0.0f }, pj_py[8] = { 0.0f }, pj_pz[8] = { 0.0f };
    f32 pj_vx[8] = { 0.0f }, pj_vy[8] = { 0.0f }, pj_vz[8] = { 0.0f };

    __m256 _mass, _viscosity;

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
    __m256 _dx, _dy, _dz;
    __m256 _dx2, _dy2, _dz2, _dx2dy2;
    __m256 _dist2, _dist;
    __m256 _dir_x, _dir_y, _dir_z;

    // Kernels
    __m256 _spiky_grad, _spiky_laplacian, _support_radius, _support_radius2;
    __m256 _dr, _dr2, _dr3;

    // Forces
    __m256 _acceleration_x, _acceleration_y, _acceleration_z;
    __m256 _gravity_force_x, _gravity_force_y, _gravity_force_z;

    // Pressure
    __m256 _pressure_force_x, _pressure_force_y, _pressure_force_z;
    __m256 _pc, _pb, _pe, _pa, _pd, _pf;
    
    // Viscosity
    __m256 _viscous_comp_x, _viscous_comp_y, _viscous_comp_z, _viscous_force_x, _viscous_force_y, _viscous_force_z;
    __m256 _dvx, _dvy, _dvz, _va, _vb, _vc;

    // Surface Tension
    __m256 _curvature, _surface_tension_x, _surface_tension_y, _surface_tension_z;
    __m256 _normal_x, _normal_y, _normal_z;
    __m256 _norm_mag, _norm_x, _norm_y, _norm_z;

    // Constants
    __m256 _zero, _two;
    __m256 _mask1, _mask2;
    
    __m256 _min, _max;

    _mass = _mm256_set1_ps(settings.mass);
    _viscosity = _mm256_set1_ps(settings.viscosity);

    _support_radius    = _mm256_set1_ps(settings.support_radius);
    _support_radius2   = _mm256_set1_ps(settings.support_radius2);
    _spiky_grad        = _mm256_set1_ps(settings.spiky_grad);
    _spiky_laplacian   = _mm256_set1_ps(settings.viscosity_laplacian);

    _gravity_x = _mm256_set1_ps(settings.gravity.x);
    _gravity_y = _mm256_set1_ps(settings.gravity.y);
    _gravity_z = _mm256_set1_ps(settings.gravity.z);

    _zero = _mm256_set1_ps(0.0f);
    _two  = _mm256_set1_ps(2.0f);

    _min = _mm256_set1_ps(std::numeric_limits<f32>::min());
    _max = _mm256_set1_ps(std::numeric_limits<f32>::max());
    auto clamp = [](__m256 a, __m256 min, __m256 max) { return _mm256_max_ps(min, _mm256_min_ps(a, max)); };

    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        Particle& pi = fluid_particles[i];

        // Set pi data
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
        _pf               = _mm256_setzero_ps();
        _pressure_force_x = _mm256_setzero_ps();
        _pressure_force_y = _mm256_setzero_ps();
        _pressure_force_z = _mm256_setzero_ps();

        // Reset viscous
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

        vi3 cell_origin = cell(pi, settings.support_radius);
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

                        // Distance squared
                        _dx     = _mm256_sub_ps(_pj_px, _pi_px);
                        _dy     = _mm256_sub_ps(_pj_py, _pi_py);
                        _dz     = _mm256_sub_ps(_pj_pz, _pi_pz);
                        _dx2    = _mm256_mul_ps(_dx, _dx);
                        _dy2    = _mm256_mul_ps(_dy, _dy);
                        _dz2    = _mm256_mul_ps(_dz, _dz);
                        _dx2dy2 = _mm256_add_ps(_dx2, _dy2);
                        _dist2  = _mm256_add_ps(_dx2dy2, _dz2);

                        // if dist2 < support_radius
                        _mask1 = _mm256_cmp_ps(_dist2, _support_radius2, _CMP_LT_OS);

                        // Distance
                        //dist = sqrt(dist2);
                        _dist = _mm256_sqrt_ps(_dist2);

                        // TODO: REINPLEMENT PRESSURE FORCE
                        // Pressure force
                        // pressure_force += pj.mass * diff * ((pi.pressure + pj.pressure) / (2.0f * pj.density)) * settings.spiky_grad * std::pow(settings.support_radius - dist, 3.0f) / dist;

                        // (support_radius - distance)^3
                        _dr  = _mm256_sub_ps(_support_radius, _dist);
                        _dr2 = _mm256_mul_ps(_dr, _dr);
                        _dr3 = _mm256_mul_ps(_dr2, _dr);

                        // pa = pi.pressure + pj.pressure
                        _pa = _mm256_add_ps(_pi_pressure, _pj_pressure);
                        // pb = 2 * pj.density
                        _pb = _mm256_mul_ps(_two, _pj_density);
                        // pc = ((pi.pressure + pj.pressure) / (2 * pj.density))
                        _pc = _mm256_div_ps(_pa, _pb);
                        // pd = pj.mass * ((pi.pressure + pj.pressure)  / (2 * pj.density)) 
                        _pd = _mm256_mul_ps(_mass, _pc);
                        // pe = pj.mass * ((pi.pressure + pj.pressure)  / (2 * pj.density)) * spiky_grad
                        _pe = _mm256_mul_ps(_pd, _spiky_grad);
                        // pf = pj.mass * ((pi.pressure + pj.pressure)  / (2 * pj.density)) * spiky_grad * (support_radius - distance)^3
                        _pf = _mm256_mul_ps(_pe, _dr3);
                        // if dist2 < radius2
                        _pf = _mm256_blendv_ps(_zero, _pf, _mask1);

                        // direction = glm::normalize(difference);
                        _dir_x = _mm256_div_ps(_dx, _dist);
                        _dir_y = _mm256_div_ps(_dy, _dist);
                        _dir_z = _mm256_div_ps(_dz, _dist);

                        // pressure_force += -direction * pj.mass * ((pi.pressure + pj.pressure) / (2.0f * pj.density)) * spiky_grad * (support_radius - distance)^3
                        _pressure_force_x = _mm256_fmadd_ps(_pf, _dir_x, _pressure_force_x);
                        _pressure_force_y = _mm256_fmadd_ps(_pf, _dir_y, _pressure_force_y);
                        _pressure_force_z = _mm256_fmadd_ps(_pf, _dir_z, _pressure_force_z);

                        // TODO: REINPLEMENT VISCOUS FORCE
                        // Viscous force
                        // viscous_force += pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density) * spiky_laplacian * (support_radius - distance)

                        // (pj.velocity - pi.velocity)
                        _dvx = _mm256_sub_ps(_pj_vx, _pi_vx);
                        _dvy = _mm256_sub_ps(_pj_vy, _pi_vy);
                        _dvz = _mm256_sub_ps(_pj_vz, _pi_vz);

                        // ((pj.velocity - pi.velocity) / pj.density)
                        _va = _mm256_div_ps(_dvx, _pj_density);
                        _vb = _mm256_div_ps(_dvy, _pj_density);
                        _vc = _mm256_div_ps(_dvz, _pj_density);

                        // viscosity * ((pj.velocity - pi.velocity) / pj.density)
                        _viscous_comp_x = _mm256_mul_ps(_va, _viscosity);
                        _viscous_comp_y = _mm256_mul_ps(_vb, _viscosity);
                        _viscous_comp_z = _mm256_mul_ps(_vc, _viscosity);

                        // pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density)
                        _viscous_comp_x = _mm256_mul_ps(_viscous_comp_x, _mass);
                        _viscous_comp_y = _mm256_mul_ps(_viscous_comp_y, _mass);
                        _viscous_comp_z = _mm256_mul_ps(_viscous_comp_z, _mass);

                        // pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density) * spiky_laplacian
                        _viscous_comp_x = _mm256_mul_ps(_viscous_comp_x, _spiky_laplacian);
                        _viscous_comp_y = _mm256_mul_ps(_viscous_comp_y, _spiky_laplacian);
                        _viscous_comp_z = _mm256_mul_ps(_viscous_comp_z, _spiky_laplacian);

                        // pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density) * spiky_laplacian * (support_radius - distance)
                        _viscous_comp_x = _mm256_mul_ps(_viscous_comp_x, _dr);
                        _viscous_comp_y = _mm256_mul_ps(_viscous_comp_y, _dr);
                        _viscous_comp_z = _mm256_mul_ps(_viscous_comp_z, _dr);

                        // Apply the mask to the viscous force components
                        _viscous_comp_x = _mm256_blendv_ps(_zero, _viscous_comp_x, _mask1);
                        _viscous_comp_y = _mm256_blendv_ps(_zero, _viscous_comp_y, _mask1);
                        _viscous_comp_z = _mm256_blendv_ps(_zero, _viscous_comp_z, _mask1);

                        // viscous_force += pj.mass * viscosity * ((pj.velocity - pi.velocity) / pj.density) * spiky_laplacian * (support_radius - distance);
                        _viscous_force_x = _mm256_add_ps(_viscous_force_x, _viscous_comp_x);
                        _viscous_force_y = _mm256_add_ps(_viscous_force_y, _viscous_comp_y);
                        _viscous_force_z = _mm256_add_ps(_viscous_force_z, _viscous_comp_z);

                        // TODO: REINPLEMENT SURFACE TENSION FORCE
                        // Normal
                        _normal_x = _mm256_fmadd_ps(_pf, _dir_x, _normal_x);
                        _normal_y = _mm256_fmadd_ps(_pf, _dir_y, _normal_y);
                        _normal_z = _mm256_fmadd_ps(_pf, _dir_z, _normal_z);
                    }
                }
            }
        }
        
        // TODO: REINPLEMENT SURFACE TENSION FORCE
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

        // Normalize the normal vector
        _norm_x = _mm256_div_ps(_normal_x, _norm_mag);
        _norm_y = _mm256_div_ps(_normal_y, _norm_mag);
        _norm_z = _mm256_div_ps(_normal_z, _norm_mag);

        // Surface tension curvature
        _curvature = _mm256_sub_ps(_zero, _norm_mag); // Negative curvature indicates a concave surface
        _mask2 = _mm256_cmp_ps(_curvature, _zero, _CMP_GT_OS);

        // Surface tension forces
        _surface_tension_x = _mm256_blendv_ps(_surface_tension_x, _mm256_mul_ps(_curvature, _norm_x), _mask2);
        _surface_tension_y = _mm256_blendv_ps(_surface_tension_y, _mm256_mul_ps(_curvature, _norm_y), _mask2);
        _surface_tension_z = _mm256_blendv_ps(_surface_tension_z, _mm256_mul_ps(_curvature, _norm_z), _mask2);

        // Update particle acceleration
        // pi.acceleration += (pressure_force + viscous_force + surface_tension_force) / pi.density + gravity_force;

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
        f32 acceleration_x[8] = { 0.0f }, acceleration_y[8] = { 0.0f }, acceleration_z[8] = { 0.0f };
        _mm256_storeu_ps(acceleration_x, _acceleration_x);
        _mm256_storeu_ps(acceleration_y, _acceleration_y);
        _mm256_storeu_ps(acceleration_z, _acceleration_z);

        f32 acc_x_sum = 0.0f, acc_y_sum = 0.0f, acc_z_sum = 0.0f;
        for (s32 j = 0; j < 8; j++)
        {
            if (std::isfinite(acceleration_x[j])) acc_x_sum += acceleration_x[j];
            if (std::isfinite(acceleration_y[j])) acc_y_sum += acceleration_y[j];
            if (std::isfinite(acceleration_z[j])) acc_z_sum += acceleration_z[j];
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
        vi3 cell_origin = cell(pi, settings.support_radius);
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

                        vf3 dp  = pj.position - pi.position;
                        f32 dist2 = glm::length2(dp);
                        f32 dist  = std::sqrt(dist2);

                        if (dist < settings.radius)
                        {
                            if (dist < 1e-6f)
                                dist = 1e-6f;

                            // normalize
                            vf3 dir = dp / dist;
                            f32 overlap = settings.radius - dist;

                            pi.position -= dir * overlap * 0.5f;
                            pj.position += dir * overlap * 0.5f;

                            // TODO: transform into optimal form
                            
                            //vf3 surface_normal = -diff / dist;
                            //pi.velocity += glm::dot(pj.velocity, surface_normal) * 0.005f * 2.0f;
                            //pj.velocity += glm::dot(pi.velocity, surface_normal) * 0.005f * 2.0f;

                            vf3 relative_velocity = pi.velocity - pj.velocity;
                            f32 velocity_along_normal = glm::dot(relative_velocity, dir);
                            f32 combined_mass = pi.mass + pj.mass;
                            f32 impulse = (2.0f * velocity_along_normal) / combined_mass;
                            pi.velocity -= impulse * (pj.mass / combined_mass) * dir * 0.005f;
                            pj.velocity += impulse * (pi.mass / combined_mass) * dir * 0.005f;
                        }
                    }
                }
            }
        }

        // Boundary Collision
        if (pi.position.x - settings.radius < settings.boundary_min.x)
        {
            pi.position.x = settings.boundary_min.x + settings.radius;
            pi.velocity.x *= settings.boundary_damping;
        }
        if (pi.position.x + settings.radius > settings.boundary_max.x)
        {
            pi.position.x = settings.boundary_max.x - settings.radius;
            pi.velocity.x *= settings.boundary_damping;
        }

        if (pi.position.y - settings.radius < settings.boundary_min.y)
        {
            pi.position.y = settings.boundary_min.y + settings.radius;
            pi.velocity.y *= settings.boundary_damping;
        }
        if (pi.position.y + settings.radius > settings.boundary_max.y)
        {
            pi.position.y = settings.boundary_max.y - settings.radius;
            pi.velocity.y *= settings.boundary_damping;
        }

        if (pi.position.z - settings.radius < settings.boundary_min.z)
        {
            pi.position.z = settings.boundary_min.z + settings.radius;
            pi.velocity.z *= settings.boundary_damping;
        }
        if (pi.position.z + settings.radius > settings.boundary_max.z)
        {
            pi.position.z = settings.boundary_max.z - settings.radius;
            pi.velocity.z *= settings.boundary_damping;
        }
    }
}

std::vector<Particle>& SPH::GetParticles() { return fluid_particles; }
SPHStatistics& SPH::GetStats() { return stats; }

void SPH::ComputeStatistics()
{
    stats.densities.clear();
    stats.densities.resize(fluid_particles.size());
    for (s32 i = 0; i < fluid_particles.size(); i++)
    {
        const Particle& p  = fluid_particles[i];
        stats.densities[i] = p.density;
    }
}

// Spatial Hash Map
void SPH::CreateNeighborTable()
{
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

vi3 SPH::cell(const Particle& p, f32 h) { return { p.position.x / h, p.position.y / h, p.position.z / h }; }
u32 SPH::hash(const vi3& cell) { return (static_cast<u32>(cell.x * 92837111) ^ static_cast<u32>(cell.y * 689287499) ^ static_cast<u32>(cell.z * 283823481)) % TABLE_SIZE; }
