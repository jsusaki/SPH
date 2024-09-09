#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 color;

out vec3 vs_position;
out vec3 vs_normal;
out vec4 vs_color;

uniform mat4 model;
uniform mat4 proj_view;

uniform float color_factor;

/*
// viridis: https://waldyrious.net/viridis-palette-generator/
// https://bids.github.io/colormap/
const int N_COLORS = 16;
const vec3 color_palette[N_COLORS] = vec3[N_COLORS](
    vec3(0.2667, 0.0039, 0.3294),   //  68,   1, 84
    vec3(0.2824, 0.1020, 0.4235),   //  72,  26, 108
    vec3(0.2784, 0.1843, 0.4902),   //  71,  47, 125
    vec3(0.2549, 0.2667, 0.5294),   //  65,  68, 135
    vec3(0.2235, 0.3373, 0.5490),   //  57,  86, 140
    vec3(0.1922, 0.4078, 0.5569),   //  49, 104, 142
    vec3(0.1647, 0.4706, 0.5569),   //  42, 120, 142
    vec3(0.1373, 0.5333, 0.5569),   //  35, 136, 142
    vec3(0.1216, 0.5961, 0.5451),   //  31, 152, 139
    vec3(0.1333, 0.6588, 0.5176),   //  34, 168, 132
    vec3(0.2078, 0.7176, 0.4745),   //  53, 183, 121
    vec3(0.3294, 0.7725, 0.4078),   //  84, 197, 104
    vec3(0.4784, 0.8196, 0.3176),   // 122, 209,  81
    vec3(0.6471, 0.8588, 0.2118),   // 165, 219,  54
    vec3(0.8235, 0.8863, 0.1059),   // 210, 226,  27
    vec3(0.9922, 0.9059, 0.1451)    // 253, 231,  37
);
*/

const int N_COLORS = 32;
const vec3 color_palette[N_COLORS] = vec3[N_COLORS](
    vec3(0.2667, 0.0039, 0.3294),   // 68,  1, 84
    vec3(0.2784, 0.0510, 0.3765),   // 71, 13, 96
    vec3(0.2824, 0.0941, 0.4157),   // 72, 24, 106
    vec3(0.2824, 0.1373, 0.4549),   // 72, 35, 116
    vec3(0.2784, 0.1804, 0.4863),   // 71, 46, 124
    vec3(0.2706, 0.2196, 0.5098),   // 69, 56, 130
    vec3(0.2588, 0.2549, 0.5255),   // 66, 65, 134
    vec3(0.2431, 0.2902, 0.5373),   // 62, 74, 137
    vec3(0.2275, 0.3294, 0.5490),   // 58, 84, 140
    vec3(0.2118, 0.3647, 0.5529),   // 54, 93, 141
    vec3(0.1961, 0.3961, 0.5569),   // 50, 101, 142
    vec3(0.1804, 0.4275, 0.5569),   // 46, 109, 142
    vec3(0.1686, 0.4588, 0.5569),   // 43, 117, 142
    vec3(0.1569, 0.4902, 0.5569),   // 40, 125, 142
    vec3(0.1451, 0.5176, 0.5569),   // 37, 132, 142
    vec3(0.1333, 0.5490, 0.5529),   // 34, 140, 141
    vec3(0.1216, 0.5804, 0.5490),   // 31, 148, 140
    vec3(0.1176, 0.6118, 0.5373),   // 30, 156, 137
    vec3(0.1255, 0.6392, 0.5255),   // 32, 163, 134
    vec3(0.1451, 0.6706, 0.5098),   // 37, 171, 130
    vec3(0.1804, 0.7020, 0.4863),   // 46, 179, 124
    vec3(0.2275, 0.7294, 0.4627),   // 58, 186, 118
    vec3(0.2824, 0.7569, 0.4314),   // 72, 193, 110
    vec3(0.3451, 0.7804, 0.3961),   // 88, 199, 101
    vec3(0.4235, 0.8039, 0.3529),   // 108, 205, 90
    vec3(0.4980, 0.8275, 0.3059),   // 127, 211, 78
    vec3(0.5765, 0.8431, 0.2549),   // 147, 215, 65
    vec3(0.6588, 0.8588, 0.2039),   // 168, 219, 52
    vec3(0.7529, 0.8745, 0.1451),   // 192, 223, 37
    vec3(0.8353, 0.8863, 0.1020),   // 213, 226, 26
    vec3(0.9176, 0.8980, 0.1020),   // 234, 229, 26
    vec3(0.9922, 0.9059, 0.1451)    // 253, 231, 37
);

// TODO: oklab viridis: https://bottosson.github.io/posts/oklab/

vec3 interpolate_color(float t) 
{
    float segment_length = 1.0 / (N_COLORS - 1);
    float segment = t / segment_length;
    int idx       = int(segment);
    float frac    = segment - float(idx);

    if (idx >= N_COLORS - 1)
    {
        idx = N_COLORS - 2;
        frac = 1.0;
    }

    vec3 start_color = color_palette[idx];
    vec3 end_color   = color_palette[idx + 1];

    return mix(start_color, end_color, frac);
}

void main()
{
	vs_position = vec3(model * vec4(position, 1.0));
	vs_normal   = mat3(transpose(inverse(model))) * normal;
	gl_Position = proj_view * vec4(vs_position, 1.0f);

	float normalized_color_factor = clamp(color_factor, 0.0, 1.0);
    vec3 gradient_color = interpolate_color(normalized_color_factor);
	vs_color = vec4(mix(color.rgb, gradient_color, 1.0), color.a);	
}