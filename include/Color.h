/*
	Color Theory and Color Spaces
		viridis palette generator: https://waldyrious.net/viridis-palette-generator/
		viridis colormap: https://bids.github.io/colormap/
		colorspacious: https://github.com/njsmith/colorspacious
		colormap: https://github.com/bids/colormap
		viscm: https://github.com/matplotlib/viscm
		oklab, A perceptual color space for image processing: https://bottosson.github.io/posts/oklab/
		Perceptually uniform color spaces: https://programmingdesignsystems.com/color/perceptually-uniform-color-spaces/
		An interactive review of Oklab: https://raphlinus.github.io/color/2021/01/18/oklab-critique.html
		Online Colorblind simulator: https://daltonlens.org/colorblindness-simulator
		coloria: https://github.com/coloria-dev/coloria

		// Oklab: https://bottosson.github.io/posts/oklab/
		// ok_color: https://bottosson.github.io/misc/ok_color.h
		struct Lab { float L; float a; float b; };
		struct RGB { float r; float g; float b; };

		Lab linear_srgb_to_oklab(RGB c)
		{
			float l = 0.4122214708f * c.r + 0.5363325363f * c.g + 0.0514459929f * c.b;
			float m = 0.2119034982f * c.r + 0.6806995451f * c.g + 0.1073969566f * c.b;
			float s = 0.0883024619f * c.r + 0.2817188376f * c.g + 0.6299787005f * c.b;

			float l_ = cbrtf(l);
			float m_ = cbrtf(m);
			float s_ = cbrtf(s);

			return {
				0.2104542553f * l_ + 0.7936177850f * m_ - 0.0040720468f * s_,
				1.9779984951f * l_ - 2.4285922050f * m_ + 0.4505937099f * s_,
				0.0259040371f * l_ + 0.7827717662f * m_ - 0.8086757660f * s_,
			};
		}

		RGB oklab_to_linear_srgb(Lab c)
		{
			float l_ = c.L + 0.3963377774f * c.a + 0.2158037573f * c.b;
			float m_ = c.L - 0.1055613458f * c.a - 0.0638541728f * c.b;
			float s_ = c.L - 0.0894841775f * c.a - 1.2914855480f * c.b;

			float l = l_ * l_ * l_;
			float m = m_ * m_ * m_;
			float s = s_ * s_ * s_;

			return {
				+4.0767416621f * l - 3.3077115913f * m + 0.2309699292f * s,
				-1.2684380046f * l + 2.6097574011f * m - 0.3413193965f * s,
				-0.0041960863f * l - 0.7034186147f * m + 1.7076147010f * s,
			};
		}
*/
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

static const std::vector<ucolor> palette = {
	//{ 255, 0, 0, 255 },   // Red
	//{ 0, 255, 0, 255 },   // Green
	//{ 0, 0, 255, 255 },   // Blue
	//{ 255, 255, 0, 255 }, // Yellow
	//{ 255, 0, 255, 255 }, // Magenta
	//{ 0, 255, 255, 255 }  // Cyan
	{ 0, 102, 178, 128 }, // water
};

// Viridis palette
static const std::vector<ucolor> viridis_16 = {
	{ 68,   1, 84},
	{ 72,  26, 108},
	{ 71,  47, 125},
	{ 65,  68, 135},
	{ 57,  86, 140},
	{ 49, 104, 142},
	{ 42, 120, 142},
	{ 35, 136, 142},
	{ 31, 152, 139},
	{ 34, 168, 132},
	{ 53, 183, 121},
	{ 84, 197, 104},
	{122, 209,  81},
	{165, 219,  54},
	{210, 226,  27},
	{253, 231,  37},
};

static const std::vector<ucolor> viridis_32 = {
	{ 68,  1, 84 },
	{ 71, 13, 96 },
	{ 72, 24, 106 },
	{ 72, 35, 116 },
	{ 71, 46, 124 },
	{ 69, 56, 130 },
	{ 66, 65, 134 },
	{ 62, 74, 137 },
	{ 58, 84, 140 },
	{ 54, 93, 141 },
	{ 50, 101, 142 },
	{ 46, 109, 142 },
	{ 43, 117, 142 },
	{ 40, 125, 142 },
	{ 37, 132, 142 },
	{ 34, 140, 141 },
	{ 31, 148, 140 },
	{ 30, 156, 137 },
	{ 32, 163, 134 },
	{ 37, 171, 130 },
	{ 46, 179, 124 },
	{ 58, 186, 118 },
	{ 72, 193, 110 },
	{ 88, 199, 101 },
	{ 108, 205, 90 },
	{ 127, 211, 78 },
	{ 147, 215, 65 },
	{ 168, 219, 52 },
	{ 192, 223, 37 },
	{ 213, 226, 26 },
	{ 234, 229, 26 },
	{ 253, 231, 37 },
};