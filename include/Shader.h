#pragma once

#include <string>
#include <print>
#include <fstream>
#include <vector>
#include <iostream>
#include <unordered_map>

#include <glad/glad.h>

#include "Common.h"


enum ShaderType
{
	VERTEX   = GL_VERTEX_SHADER,
	FRAGMENT = GL_FRAGMENT_SHADER,
	GEOMETRY = GL_GEOMETRY_SHADER,
	COMPUTE  = GL_COMPUTE_SHADER,
};

class Shader
{
public:
	Shader();
	Shader(const std::string& vertex_path, const std::string& fragment_path);
	~Shader();

public:
	void Use();
	void Unuse();
	u32 GetID();

private:
	std::string load_from_file(const std::string& filepath);
	u32  compile(ShaderType type, std::string& source);
	void attach(u32& shader);
	void link();
	void detach(u32& shader);
	void destroy(u32& shader);

public:
	u32  GetAttribute(const std::string& name) const;
	u32  GetUniform(const std::string& name)  const;
	void SetUniform(const std::string& name, const s32& val);
	void SetUniform(const std::string& name, f32* val, s32 count);
	void SetUniform(const std::string& name, s32* val, s32 count);
	void SetUniform(const std::string& name, const f64& val);
	void SetUniform(const std::string& name, const f32& val);
	void SetUniform(const std::string& name, const vf2& vector);
	void SetUniform(const std::string& name, const vf3& vector);
	void SetUniform(const std::string& name, const vf4& vector);
	void SetUniform(const std::string& name, const mf4x4& matrix);

private:
	u32 id;
	mutable std::unordered_map<std::string, u32> uniform_locations;
};