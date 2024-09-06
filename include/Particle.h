#pragma once

#include "Common.h"

struct Particle
{
    vf3 position;
    vf3 velocity;
    vf3 acceleration;
    f32 mass;
    f32 density;
    f32 pressure;
    f32 viscosity;
    f32 radius;
    f32 support_radius;
    u32 hash;
};