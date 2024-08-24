#pragma once

// C-Standard types
#include <cstdint>
// Numeric types
using u8  = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using s8  = int8_t;
using s16 = int16_t;
using s32 = int32_t;
using s64 = int64_t;
using f32 = float;
using f64 = double;

// GLM types
#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>
// Vector types
using vu2 = glm::uvec2;
using vi2 = glm::ivec2;
using vf2 = glm::fvec2;
using vd2 = glm::dvec2;
using vu3 = glm::uvec3;
using vi3 = glm::ivec3;
using vf3 = glm::fvec3;
using vd3 = glm::dvec3;
using vu4 = glm::uvec4;
using vi4 = glm::ivec4;
using vf4 = glm::fvec4;
using vd4 = glm::dvec4;
// Matrix types
using mi3x3 = glm::imat3x3;
using mf3x3 = glm::fmat3x3;
using md3x3 = glm::dmat3x3;
using mi4x4 = glm::imat4x4;
using mf4x4 = glm::fmat4x4;
using md4x4 = glm::dmat4x4;
// Quaternion types
using qf32 = glm::f32quat;
using qf64 = glm::f64quat;