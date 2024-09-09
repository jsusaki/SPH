#version 330 core

in vec3 vs_position;
in vec3 vs_normal;
in vec4 vs_color;

out vec4 fs_color;

void main()
{ 
	vec3 norm = normalize(vs_normal);
	fs_color  = vec4(vs_color);
}