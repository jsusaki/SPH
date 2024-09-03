#pragma once

#include "Common.h"

const u32 SCREEN_WIDTH = 1280;
const u32 SCREEN_HEIGHT = 720;

static const f32 PI = 3.14159265358979323846f;

namespace config
{
    // Simulation
    static const s32 NUM_PARTICLES = 5000;
    static const f32 REST_DENSITY = 997.0f;
    static const f32 GAS_CONSTANT = 1.0f;    // Gas Constant, or stiffness control the particle spread, make 5 to implode
    static const f32 VISCOSITY = 0.002f;  // Kinematic Viscosity
    static const f32 SURFACE_TENSION_CONSTANT = 0.0000001f; // 0.0728f
    static const vf3 GRAVITY = { 0.0f, -9.80665f, 0.0f };
    static const f32 MASS = 0.01f;
    static const f32 DT = 0.0001f;
    // Smoothing Kernel
    static const f32 SMOOTHING_RADIUS = 0.04f;
    static const f32 SMOOTHING_RADIUS_2 = std::pow(SMOOTHING_RADIUS, 2.0f);
    static const f32 POLY6 = 315.0f / (64.0f * PI * std::pow(SMOOTHING_RADIUS, 9.0f));   // Density
    static const f32 SPIKY_GRAD = -45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f));           // Pressure
    static const f32 SPIKY_LAPLACIAN = 45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f));           // Viscosity
    // Boundary
    static const f32 BOUNDARY_EPSILON = 0.0000001f;
    static const f32 BOUNDARY_DAMPING = -0.4f;
    static const vf3 BOUNDARY = { 0.6f, 0.6f, 0.6f };
    static const vf3 BOUNDARY_MIN = -BOUNDARY;
    static const vf3 BOUNDARY_MAX = BOUNDARY;
    // Hash map
    const u32 TABLE_SIZE = 262144; //NUM_PARTICLES * 2; //262144;
    const u32 NO_PARTICLE = 0xFFFFFFFF;
}

using namespace config;