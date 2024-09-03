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
    void ResetSimulation(const SPHSettings& settings);
    void RestartSimulation();

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
    vf2 mouse_pos_transformed = {};
    vf2 prev_mouse_pos_transformed = {};

private: // Simulation variables
    // Camera
    ArcballCamera camera;
    vf2 transform_mouse(vf2 p) { return vf2(p.x * 2.0f / SCREEN_WIDTH - 1.0f, 1.0f - 2.0f * p.y / SCREEN_HEIGHT); }

    // Graphics
    std::shared_ptr<Shader> shader;
    std::vector<model> particle_models;
    std::vector<model> boundary_models;

    // Smoothed Particle Hydrodynamics
    SPH sph;
    SPHSettings default_settings;
    SPHSettings current_settings;

    // SPH Simulation State
    enum class State : u8
    {
        PAUSE,     // Pause the simulation
        SIMULATE,  // Run the simulation
        STEP,      // Step the simulation
        RESET,     // Reset the simulation default state
        APPLY,     // Apply new settings
        RESTART,   // Restart the current simulation
    };
    State current_state = State::PAUSE;
    State next_state    = State::PAUSE;
    s32 n_steps = 0;

    // GUI
    GUI gui;
};
