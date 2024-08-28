#pragma once

#include <print>
#include <memory>
#include <chrono>
#include <algorithm>

#include "Common.h"
#include "Random.h"
#include "SPH.h"

#include "Shader.h"
#include "Camera.h"
#include "Mesh.h"
#include "Window.h"
#include "Input.h"
#include "GUI.h"

#include "../test/Test.h"

class Simulator
{
public:
    Simulator();

public: // Interface
    bool Init();
    bool Start();
    bool ShutDown();

private: // Main functions
    void Create();
    void ProcessInput();
    void Simulate(f32 dt);
    void Render();

private: // Helper functions
    void UpdateFrameTime();

private: // Simulator variables
    // Timing
    std::chrono::time_point<std::chrono::system_clock> m_t1;
    std::chrono::time_point<std::chrono::system_clock> m_t2;
    u32 m_last_fps;
    u32 m_frame_count;
    f32 m_frame_timer;
    f32 m_accumulator;
    f32 m_delta_time;
    f32 m_elapsed_time;
    f32 m_last_elapsed_time;

    // Window
    Window window;

    // Input
    vf2 mouse_pos;
    vf2 prev_mouse_pos;
    bool mouse_control = true;

private: // Simulation variables
    // Camera
    std::shared_ptr<ArcballCamera> camera;
    vf2 transform_mouse(vf2 p) { return vf2(p.x * 2.0f / SCREEN_WIDTH - 1.0f, 1.0f - 2.0f * p.y / SCREEN_HEIGHT); }

    // Graphics
    std::shared_ptr<Shader> shader;
    std::vector<model> particle_models;
    std::vector<model> boundary_models;

    // Smoothed Particle Hydrodynamics
    SPH sph;
    settings sph_settings;
    bool simulate = false;
    bool reset = false;
    bool step = false;

    // GUI
    GUI gui;

};
