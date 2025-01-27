#pragma once

#include "Common.h"
#include "SPH.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

class ArcBallCamera;

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

        // Setup ImPlot context
        ImPlot::CreateContext();
	}

    void Display(SPHSettings& s, ArcballCamera& c, SPHStatistics& stats)
    {
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        {
            ImGui::Begin("SPH Simulator");

            // General
            ImGui::SeparatorText("General");

            // Camera position
            ImGui::Text("Eye:    x=%.3f y=%.3f z=%.3f", c.eye().x, c.eye().y, c.eye().z);
            ImGui::Text("Up:     x=%.3f y=%.3f z=%.3f", c.up().x,  c.up().y,  c.up().z);
            ImGui::Text("Center: x=%.3f y=%.3f z=%.3f", c.center().x, c.center().y, c.center().z);

            // Frame Rate
            ImGui::Text("FPS: average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

            ImGui::Checkbox("Wireframe", &enable_wireframe);

            // SPH Settings
            ImGui::SeparatorText("Settings");
            static ImGuiSliderFlags slider_flags = ImGuiSliderFlags_AlwaysClamp;
            ImGui::SliderInt("Number of Particles",        &s.n_particles,                  1,         25000);
            ImGui::SliderFloat("Rest Density",             &s.rest_density,                 0.001f,    1000.0f);
            ImGui::SliderFloat("Stiffness",                &s.stiffness,                    0.001f,    10.0f);
            ImGui::SliderFloat("Kinematic Viscosity",      &s.viscosity,                    0.001f,    10.0f);
            ImGui::SliderFloat("Surface Tension Constant", &s.surface_tension_constant,     0.000001f, 1.0f);
            ImGui::SliderFloat("Mass",                     &s.mass,                         0.001f,    10.0f);
            ImGui::SliderFloat3("Gravity",                 glm::value_ptr(s.gravity),       -10.0f,    10.0f);
            ImGui::SliderFloat("Support Radius",           &s.support_radius,               0.001f,    0.1f);
            ImGui::SliderFloat("Time Step",                &s.dt,                           0.0001f,   0.5f);
            ImGui::SliderFloat("Boundary Damping",         &s.boundary_damping,             -1.0f,     1.0f);
            ImGui::SliderFloat3("Boudary Size",            glm::value_ptr(s.boundary_size), 0.1f,      1.0f);
            ImGui::SliderFloat("Particle Spacing",         &s.particle_spacing,             0.0f,      1.0f);
            ImGui::SliderFloat3("Cube Offset",             glm::value_ptr(s.cube_offset),   0.0f,      1.0f);

            apply_pressed = false; 
            if (ImGui::Button("Apply")) apply_pressed = true;     ImGui::SameLine();

            restart_pressed = false;
            if (ImGui::Button("Restart")) restart_pressed = true; ImGui::SameLine();

            reset_pressed = false;
            if (ImGui::Button("Reset")) reset_pressed = true;

            ImGui::Checkbox("AVX2",             &enable_avx2);
            ImGui::Checkbox("Show Density Distribution", &show_density_distribution);
            ImGui::Checkbox("Show ImGui Demo",  &show_imgui_demo);
            ImGui::Checkbox("Show ImPlot Demo", &show_implot_demo);

            if (show_imgui_demo)  ImGui::ShowDemoWindow();
            if (show_implot_demo) ImPlot::ShowDemoWindow();

            // Particle Density Distribution
            if (show_density_distribution)
            {
                ImGui::Begin("SPH Density Distribution");

                static s32 xmin = 24470, xmax = 24800;
                static f32 ymin = 0.0, ymax = 0.15;
                static ImPlotHistogramFlags hist_flags = ImPlotHistogramFlags_Density;
                ImPlot::SetNextAxesLimits(xmin, xmax, ymin, ymax, ImGuiCond_Always);
                if (ImPlot::BeginPlot("Particle Density Distribution"))
                {
                    // ImPlotAxisFlags_RangeFit, ImPlotAxisFlags_AutoFit
                    ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
                    ImPlot::PlotHistogram("Density", stats.densities.data(), stats.densities.size(), stats.n_bins, 1.0, ImPlotRange(xmin, xmax), hist_flags);
                    //ImPlot::PlotHistogram("Density", stats.densities.data(), stats.densities.size(), stats.n_bins, 1.0, ImPlotRange() /*ImPlotRange(xmin, xmax)*/, hist_flags);
                }

                ImGui::SliderInt("X Min",   &xmin, 24000, 25000);
                ImGui::SliderInt("X Max",   &xmax, 24000, 25000);
                ImGui::SliderFloat("Y Min", &ymin, 0.0, 1.0);
                ImGui::SliderFloat("Y Max", &ymax, 0.0, 1.0);

                ImPlot::EndPlot();
            }

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
        ImPlot::DestroyContext();
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
    }

    bool IsWindowFocused() const { return ImGui::IsWindowFocused(ImGuiFocusedFlags_AnyWindow); }

    bool IsApplyPressed()   const { return apply_pressed; }
    bool IsResetPressed()   const { return reset_pressed; }
    bool IsRestartPressed() const { return restart_pressed; }

    bool IsAVX2Enabled()    const { return enable_avx2; }
    void ToggleAVX2()  { enable_avx2 = !enable_avx2; }

    bool IsWireframeEnabled() const { return enable_wireframe; }
    void ToggleWireframe() { enable_wireframe = !enable_wireframe; }

private:
    ImGuiIO io;
    bool apply_pressed = false;
    bool reset_pressed = false;
    bool restart_pressed = false;
    bool enable_avx2 = false;
    bool enable_wireframe = false;

    bool show_density_distribution = false;

    // demo reference
    bool show_imgui_demo = false;
    bool show_implot_demo = false;

};