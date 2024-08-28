#pragma once

#include <algorithm>
#include <vector>

#include "Common.h"

// float 32-bit Color 
// TODO: vf4
struct fcolor
{
	f32 r, g, b, a;
};

// uint 32-bit Color 
// TODO: vu4
struct ucolor
{
	union { u32 p; struct { u8 r; u8 g; u8 b; u8 a; }; };
	ucolor(u8 r, u8 g, u8 b, u8 a = 0xFF) : r(r), g(g), b(b), a(a) {}
};

static fcolor to_float(const ucolor& c)
{
	return fcolor(c.r / 255.0f, c.g / 255.0f, c.b / 255.0f, c.a / 255.0f);
}

static ucolor to_uint(const fcolor& c)
{
	u8 r = static_cast<u8>(std::clamp(c.r, 0.0f, 1.0f) * 255.0f);
	u8 g = static_cast<u8>(std::clamp(c.g, 0.0f, 1.0f) * 255.0f);
	u8 b = static_cast<u8>(std::clamp(c.b, 0.0f, 1.0f) * 255.0f);
	u8 a = static_cast<u8>(std::clamp(c.a, 0.0f, 1.0f) * 255.0f);
	return ucolor(r, g, b, a);
}

static const std::vector<ucolor> colors = {
	//{ 255, 0, 0, 255 },   // Red
	//{ 0, 255, 0, 255 },   // Green
	//{ 0, 0, 255, 255 },   // Blue
	//{ 255, 255, 0, 255 }, // Yellow
	//{ 255, 0, 255, 255 }, // Magenta
	//{ 0, 255, 255, 255 }  // Cyan

	{ 0, 102, 178, 128 }, // water
};