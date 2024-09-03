#pragma once

#include "Common.h"
#include "SPH.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

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
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
        // Setup Dear ImGui style
        ImGui::StyleColorsDark();
        // Setup Platform/Renderer backends
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 130");
        ImGui::SetNextWindowFocus();
	}

    void DisplaySettings(SPHSettings& s)
    {
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            ImGui::Begin("SPH Settings");

            ImGui::Text("Change simulation parameters. Press Reset to apply changes.");

            ImGui::SliderInt("Number of Particles",        &s.n_particles,                  1000,   25000);
            ImGui::SliderFloat("Rest Density",             &s.rest_density,                 0.001f, 1000.0f);
            ImGui::SliderFloat("Gas Constant",             &s.gas_constant,                 0.001f, 10.0f);
            ImGui::SliderFloat("Kinematic Viscosity",      &s.viscosity,                    0.001f, 10.0f);
            ImGui::SliderFloat("Surface Tension Constant", &s.surface_tension_constant,     0.0000001f, 0.1f);
            ImGui::SliderFloat("Mass",                     &s.mass,                         0.001f, 10.0f);
            ImGui::SliderFloat3("Gravity",                 glm::value_ptr(s.gravity),       -10.0f, 10.0f);
            ImGui::SliderFloat("Smoothing Radius",         &s.smoothing_radius,             0.001f, 0.1f);
            ImGui::SliderFloat("Boundary Damping",         &s.boundary_damping,             -1.0f,   1.0f);
            ImGui::SliderFloat3("Boudary Size",            glm::value_ptr(s.boundary_size), 0.1f,   1.0f);

            if (ImGui::Button("Apply"))   apply_pressed = true;
            if (ImGui::Button("Reset"))   reset_pressed = true;
            if (ImGui::Button("Restart")) restart_pressed = true;

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }
    }

    void Render()
    {
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    void Shutdown()
    {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
    }

    bool IsWindowFocused() const { return ImGui::IsWindowFocused(ImGuiFocusedFlags_AnyWindow); }

    bool IsApplyPressed() const { return apply_pressed; }
    void ClearApplyRequest() { apply_pressed = false; }
    bool IsResetPressed() const { return reset_pressed; }
    void ClearResetRequest() { reset_pressed = false; }
    bool IsRestartPressed() const { return restart_pressed; }
    void ClearRestartRequest() { restart_pressed = false; }

private:
    ImGuiIO io;
    bool apply_pressed = false;
    bool reset_pressed = false;
    bool restart_pressed = false;
};