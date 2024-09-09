#include "Shader.h"

Shader::Shader()
{

}

Shader::Shader(const std::string& vertex_path, const std::string& fragment_path)
{
    id = glCreateProgram();

    std::string vertex_source   = load_from_file(vertex_path);
    std::string fragment_source = load_from_file(fragment_path);

    u32 vertex_shader   = compile(VERTEX, vertex_source);
    u32 fragment_shader = compile(FRAGMENT, fragment_source);

    attach(vertex_shader);
    attach(fragment_shader);

    link();

    detach(vertex_shader);
    detach(fragment_shader);

    destroy(vertex_shader);
    destroy(fragment_shader);
}

Shader::~Shader()
{
    glDeleteProgram(id);
}

void Shader::Use()
{
    glUseProgram(id);
}

void Shader::Unuse()
{
    glUseProgram(0);
}

u32 Shader::GetID()
{
    return id;
}

std::string Shader::load_from_file(const std::string& filepath)
{
    std::string data;
    std::ifstream file(filepath, std::ios::in | std::ios::binary);
    if (!file.is_open())
        std::cout << "Could not open " << filepath << "\n";

    std::copy(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), std::back_inserter(data));
    file.close();

    return data;
}

u32 Shader::compile(ShaderType type, std::string& source)
{
    u32 shader = glCreateShader(type);
    const char* shader_src = source.c_str();
    glShaderSource(shader, 1, &shader_src, nullptr);
    glCompileShader(shader);

    s32 status = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    {
        s32 max_length = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &max_length);

        std::vector<char> error_log(max_length);
        glGetShaderInfoLog(shader, max_length, &max_length, &error_log[0]);

        glDeleteShader(shader);

        std::printf("%s\n", &(error_log[0]));
    }

    return shader;
}

void Shader::attach(u32& shader)
{
    glAttachShader(id, shader);
}

void Shader::link()
{
    glLinkProgram(id);

    s32 status = 0;
    glGetProgramiv(id, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        s32 max_length = 0;
        glGetProgramiv(id, GL_INFO_LOG_LENGTH, &max_length);
       
        std::vector<char> error_log(max_length);
        glGetProgramInfoLog(id, max_length, &max_length, &error_log[0]);
        
        glDeleteProgram(id);
        
        std::printf("%s\n", &(error_log[0]));
    }
}

void Shader::detach(u32& shader)
{
    glDetachShader(id, shader);
}

void Shader::destroy(u32& shader)
{
    glDeleteShader(shader);
}

u32 Shader::GetAttribute(const std::string& name) const
{
    return glGetAttribLocation(id, name.c_str());
}

u32 Shader::GetUniform(const std::string& name) const
{
    if (uniform_locations.find(name) != uniform_locations.end())
        return uniform_locations[name];

    u32 location = glGetUniformLocation(id, name.c_str());
    uniform_locations[name] = location;
    return location;
}

void Shader::SetUniform(const std::string& name, const s32& val)
{
    glUniform1i(GetUniform(name), val);
}

void Shader::SetUniform(const std::string& name, f32* val, s32 count)
{
    glUniform1fv(GetUniform(name), count, val);
}

void Shader::SetUniform(const std::string& name, s32* val, s32 count)
{
    glUniform1iv(GetUniform(name), count, val);
}

void Shader::SetUniform(const std::string& name, const f64& val)
{
    glUniform1f(GetUniform(name), val);
}

void Shader::SetUniform(const std::string& name, const f32& val)
{
    glUniform1f(GetUniform(name), val);
}

void Shader::SetUniform(const std::string& name, const vf2& vector)
{
    glUniform2f(GetUniform(name), vector.x, vector.y);
}

void Shader::SetUniform(const std::string& name, const vf3& vector)
{
    glUniform3f(GetUniform(name), vector.x, vector.y, vector.z);
}

void Shader::SetUniform(const std::string& name, const vf4& vector)
{
    glUniform4f(GetUniform(name), vector.x, vector.y, vector.z, vector.w);
}

void Shader::SetUniform(const std::string& name, const mf4x4& matrix)
{
    glUniformMatrix4fv(GetUniform(name), 1, GL_FALSE, glm::value_ptr(matrix));
}