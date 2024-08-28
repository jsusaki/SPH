#pragma once

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Common.h"

class GUI
{
public:
    GUI() {}

public:
    void Init(GLFWwindow* window)
	{
        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        io = ImGui::GetIO(); (void)io;
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
        // Setup Dear ImGui style
        ImGui::StyleColorsDark();
        // Setup Platform/Renderer backends
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 130");
	}

    void Update()
    {
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
        {
            s32 n_particles = 25000;
            f32 rest_density = 997.0f;
            f32 gas_constant = 1.0f;
            f32 viscosity = 0.002f;
            f32 surface_tension_constant = 0.0000001f;
            vf3 gravity = { 0.0f, -9.80665f, 0.0f };
            f32 mass = 0.01f;
            f32 smoothing_radius = 0.04f;
            f32 boundary_damping = -0.4f;
            vf3 boundary_size = { 0.6f, 0.6f, 0.6f };

            ImGui::Begin("SPH Settings");

            ImGui::Text("Change simulation parameters. Press Reset to apply changes.");

            ImGui::SliderInt("Number of Particles", &n_particles, 1000, 25000);
            ImGui::SliderFloat("Rest Density", &rest_density, 0.001f, 1000.0f);
            ImGui::SliderFloat("Gas Constant", &gas_constant, 0.001f, 1.0f);
            ImGui::SliderFloat("Kinematic Viscosity", &viscosity, 0.001f, 10.0f);
            ImGui::SliderFloat("Surface Tension Constant", &surface_tension_constant, 0.001f, 1.0f);
            ImGui::SliderFloat("Mass", &mass, 0.001f, 10.0f);
            ImGui::SliderFloat3("Gravity", glm::value_ptr(gravity), 0.001f, 100.0f);
            ImGui::SliderFloat("Smoothing Radius", &smoothing_radius, 0.001f, 5.0f);
            ImGui::SliderFloat("Boundary Dampening", &boundary_damping, 1.0f, 100.0f);
            ImGui::SliderFloat3("Boudary Size", glm::value_ptr(boundary_size), 0.1f, 1.0f);

            if (ImGui::Button("Reset")) {

            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }
    }

    void Render()
    {
        // Rendering
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    void Shutdown()
    {
        // Cleanup
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
    }

private:
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    ImGuiIO io;
};