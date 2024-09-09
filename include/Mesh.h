#pragma once
#include <vector>
#include <fstream>
#include <algorithm>

#include <glad/glad.h>

#include "Common.h"
#include "Color.h"

struct vertex
{
	vf3 position;
	vf3 normal;
	vf4 color;
};

struct mesh
{
	u32 vao, vbo, ibo;
	std::vector<vertex> vertices;
	std::vector<u32> indices;

	mesh() {};

	mesh(std::string filepath) 
	{ 
		load_from_file(filepath);

		// Generate vertex buffers
		glGenVertexArrays(1, &vao);
		glGenBuffers(1, &vbo);
		glGenBuffers(1, &ibo);
		// Bind vertex array object
		glBindVertexArray(vao);

		// Vertex Buffer Object
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);

		// Index Buffer Object
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(u32) * indices.size(), indices.data(), GL_STATIC_DRAW);

		// Vertex Attribute configuration
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, position));

		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, normal));

		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, color));

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}

	mesh(std::vector<vertex> v, std::vector<u32> i)
	{
		vertices = v;
		indices  = i;

		// Generate vertex buffers
		glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glGenBuffers(1, &ibo);
		// Bind vertex array object
        glBindVertexArray(vao);

        // Vertex Buffer Object
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);

		// Index Buffer Object
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(u32) * indices.size(), indices.data(), GL_STATIC_DRAW);
		
		// Vertex Attribute configuration
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, position));

		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, normal));

		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, color));

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
	}

	~mesh()
	{
		//glDeleteVertexArrays(1, &vao);
		//glDeleteBuffers(1, &vbo);
		//glDeleteBuffers(1, &ibo);
	}

	void draw(s32 mode)
	{
		glBindVertexArray(vao);
		glDrawElements(mode, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	void load_from_file(std::string filepath)
	{
		std::ifstream file(filepath);
		if (!file.is_open())
			std::cout << "File does not exist\n";

		// buffers
		std::vector<vf3> vertex_points;
		std::vector<vf3> vertex_normals;

		std::vector<vi3> vertex_indices;
		std::vector<vi2> vertex_texture_indices;
		std::vector<vi3> vertex_normal_indices;

		std::string prefix;
		while (file >> prefix)
		{
			if (prefix == "v") // Point
			{
				vf3 p;
				file >> p.x >> p.y >> p.z;
				vertex_points.push_back(p);
			}
			if (prefix == "vn") // Normal
			{
				vf3 n;
				file >> n.x >> n.y >> n.z;
				vertex_normals.push_back(n);
			}
			if (prefix == "f") // Face
			{
				// Buffers
				u8 slash;
				vi3 v_idx, v_tex_idx, v_norm_idx;
				for (s32 i = 0; i < 3; i++)
					file >> v_idx[i] >> slash >> v_tex_idx[i] >> slash >> v_norm_idx[i];

				// Obj idx starts from 1, so subtract 1
				vertex_indices.push_back(v_idx - 1);
				vertex_texture_indices.push_back(v_tex_idx - 1);
				vertex_normal_indices.push_back(v_norm_idx - 1);
			}
		}

		file.close();

		for (s32 i = 0; i < vertex_indices.size(); i++)
		{
			for (s32 j = 0; j < 3; j++)
			{
				vertex v;
				v.position = vertex_points[vertex_indices[i][j]];
				v.normal   = vertex_normals[vertex_normal_indices[i][j]];

				const fcolor& c = to_float(palette[0]);
				v.color = { c.r, c.g, c.b, c.a };

				vertices.push_back(v);
				indices.push_back(vertices.size() - 1);
			}
		}
	}
};

struct cube : public mesh
{
	cube()
	{
		vertices = {
												     // color: white
			{{-1.0f, -1.0f,-1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v0
			{{ 1.0f, -1.0f,-1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v1
			{{ 1.0f,  1.0f,-1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v2
			{{-1.0f,  1.0f,-1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v3
			{{-1.0f, -1.0f, 1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v4
			{{ 1.0f, -1.0f, 1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v5
			{{ 1.0f,  1.0f, 1.0f},{1.0f, 0.0f, 0.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v6
			{{-1.0f,  1.0f, 1.0f},{0.0f, 0.0f,-1.0f},{1.0f, 1.0f, 1.0f, 1.0f}}, //v7
		};

		indices = { 
			0,1,  // v0-v1
			1,2,  // v1-v2
			2,3,  // v2-v3
			3,0,  // v3-v0
			4,5,  // v4-v5
			5,6,  // v5-v6
			6,7,  // v6-v7
			7,4,  // v7-v4
			0,4,  // v0-v4
			1,5,  // v1-v5
			2,6,  // v2-v6
			3,7   // v3-v7
		};

		// Generate vertex buffers
		glGenVertexArrays(1, &vao);
		glGenBuffers(1, &vbo);
		glGenBuffers(1, &ibo);
		// Bind vertex array object
		glBindVertexArray(vao);

		// Vertex Buffer Object
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);

		// Index Buffer Object
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(u32) * indices.size(), indices.data(), GL_STATIC_DRAW);

		// Vertex Attribute configuration
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, position));

		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, normal));

		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, color));

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);
	}
};

struct model
{
	mf4x4 m_model;
	mesh m_mesh;
	vf3 m_position;
	vf3 m_rotation;
	vf3 m_scale;
	f32 m_angle;

	model(std::string filepath)
	{
		m_model    = mf4x4(1.0f);
		m_position = vf3(0.0f, 0.0f, 0.0f);
		m_rotation = vf3(0.0f, 0.0f, 0.0f);
		m_scale    = vf3(1.0f, 1.0f, 1.0f);
		m_angle    = 0.0f;
		m_mesh     = mesh(filepath);
	}

	model(mesh m)
	{
		m_model    = mf4x4(1.0f);
		m_position = vf3(0.0f, 0.0f, 0.0f);
		m_rotation = vf3(0.0f, 0.0f, 0.0f);
		m_scale    = vf3(1.0f, 1.0f, 1.0f);
		m_angle    = 0.0f;
		m_mesh     = m;
	}

	void draw(std::shared_ptr<Shader> shader, s32 mode = GL_TRIANGLES)
	{
		shader->Use();

		m_model = mf4x4(1.0f);
		m_model = glm::rotate(m_model, glm::radians(m_angle), { 1.0f, 0.0f, 0.0f });
		m_model = glm::rotate(m_model, glm::radians(m_angle), { 0.0f, 1.0f, 0.0f });
		m_model = glm::rotate(m_model, glm::radians(m_angle), { 0.0f, 0.0f, 1.0f });

		m_model = glm::translate(m_model, m_position);
		m_model = glm::scale(m_model, m_scale);

		shader->SetUniform("model", m_model);
		m_mesh.draw(mode);
		shader->Unuse();
	}

	void scale(vf3 s) { m_scale *= s; }
	void rotate(f32 angl, vf3 axis) { m_rotation += axis; m_angle += angl; }
	void translate(vf3 pos) { m_position = pos; }
};