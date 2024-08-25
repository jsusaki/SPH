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

#include "../test/Test.h"


class Simulator
{
public:
    Simulator(): window("SPH Simulator", SCREEN_WIDTH, SCREEN_HEIGHT) {}

public:
    bool Init()
    {
        // Create user resources
        Create();

        // Time
        m_t1 = std::chrono::system_clock::now();
        m_t2 = std::chrono::system_clock::now();
        m_last_fps = 0;
        m_frame_timer = 1.0f;
        m_frame_count = 0;
        m_accumulator = 0.0f;
        m_delta_time = 1.0f / 60.0f;
        m_elapsed_time = 0.0f;
        m_last_elapsed_time = 0.0f;

        return true;
    }

    bool Start()
    {
        while (!window.ShouldClose())
        {
            // Poll events
            window.PollEvents();

            // Handle timing
            m_t2 = std::chrono::system_clock::now();
            std::chrono::duration<f32> elapsed_time = m_t2 - m_t1;
            m_t1 = m_t2;

            // Compute elapsed time
            m_elapsed_time = elapsed_time.count();
            m_last_elapsed_time = m_elapsed_time;

            // Handle input
            ProcessInput();

            // Fixed Time Update
            m_accumulator += m_delta_time;
            while (m_accumulator >= m_delta_time)
            {
                Simulate(m_elapsed_time);
                m_accumulator -= m_delta_time;
            }

            // Rendering pipeline
            Render();

            // Update Frame Time Info
            UpdateFrameTime();
        }

        return true;
    }

    bool ShutDown()
    {
        window.Close();
        return true;
    }

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

private: // Simulation variables
    // Camera
    std::shared_ptr<ArcballCamera> camera;
    mf4x4 proj;
    mf4x4 proj_inv;
    vf2 transform_mouse(vf2 p) { return vf2(p.x * 2.0f / SCREEN_WIDTH - 1.0f, 1.0f - 2.0f * p.y / SCREEN_HEIGHT); }

    vf3 eye    = { 1.5f, 0.0f, 0.0f };
    vf3 center = { 0.0f, 0.0f, 0.0f };
    vf3 up     = { 0.0f, 1.0f, 0.0f };
    f32 fov    = 65.0f;
    f32 far    = 100.0f;
    f32 near   = 0.1f;

    // Graphics
    std::shared_ptr<Shader> shader;
    std::vector<model> particle_models;
    std::vector<model> boundary_models;

    // Smoothed Particle Hydrodynamics
    SPH sph;
    bool simulate = false;
    bool reset    = false;
    bool step     = false;

private:
    void Create()
    {
        camera   = std::make_shared<ArcballCamera>(eye, center, up);
        // TODO: incorporate proj and proj_inv into camera
        proj     = glm::perspective(glm::radians(fov), static_cast<float>(SCREEN_WIDTH) / SCREEN_HEIGHT, near, far);
        proj_inv = glm::inverse(proj);

        // Create Particle Models
        mesh sphere_mesh("./res/models/sphere.obj");
        std::vector<SPH::particle>& particles = sph.GetParticles();
        for (auto& p : particles)
        {
            model sphere_model = model(sphere_mesh);
            sphere_model.translate(p.position);
            sphere_model.scale({ 0.01f, 0.01f, 0.01f });
            particle_models.push_back(sphere_model);
        }
 
        // Create Bounding Box Model
        model cube_model = model(cube());
        cube_model.scale({ config::BOUNDARY.x, config::BOUNDARY.y, config::BOUNDARY.z });
        boundary_models.push_back(cube_model);

        // Create Shader
        shader = std::make_shared<Shader>("res/shaders/basic.vs", "res/shaders/basic.fs");
    }

    void ProcessInput()
    {
        // Key control
        if (input.IsKeyPressed(GLFW_KEY_ESCAPE))
            window.SetShouldClose();

        if (input.IsKeyPressed(GLFW_KEY_SPACE))
            simulate = !simulate;

        if (input.IsKeyHeld(GLFW_KEY_S))
            step = true;

        if (input.IsKeyPressed(GLFW_KEY_R))
            reset = true;

        // Mouse control
        mouse_pos = transform_mouse(input.mouse_pos);
        if (prev_mouse_pos != (vf2(-2.0f)))
        {
            if (input.IsButtonPressed(GLFW_MOUSE_BUTTON_LEFT))
            {
                camera->rotate(prev_mouse_pos, mouse_pos);
            }
            if (input.IsButtonPressed(GLFW_MOUSE_BUTTON_RIGHT))
            {
                // TODO: incorporate into camera
                vf2 dxy  = mouse_pos - prev_mouse_pos;
                vf4 dxy4 = proj_inv * vf4(dxy.x, dxy.y, 0.0f, 1.0f);
                camera->pan(vf2(dxy4.x, dxy4.y));
            }
        }
        prev_mouse_pos = mouse_pos;

        if (input.mouse_wheel_delta != 0)
        {
            camera->zoom(input.mouse_wheel_delta * 0.1f);
            input.mouse_wheel_delta = 0;
        }

        // Window resized event
        //proj = glm::perspective(glm::radians(65.0f), static_cast<float>(SCREEN_WIDTH) / SCREEN_HEIGHT, 0.1f, 500.0f);

        input.Update();
    }

    void Simulate(f32 dt)
    {
        if (reset)
        {
            sph.Init();
            std::vector<SPH::particle>& particles = sph.GetParticles();
            for (s32 i = 0; i < particles.size(); i++)
                particle_models[i].translate(particles[i].position);
            reset    = false;
            simulate = false;
        }

        if (simulate || step)
        {
            sph.Simulate(dt);
            std::vector<SPH::particle>& particles = sph.GetParticles();
            for (s32 i = 0; i < particles.size(); i++)
                particle_models[i].translate(particles[i].position);
            if (step) step = false;
        }
    }

    void Render()
    {
        window.Clear({ 25, 25, 25, 255 });

        // Rendering 
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_MULTISAMPLE);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_POLYGON_SMOOTH);
        glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

        // Draw
        shader->Use();
        // TODO: camera get projection
        shader->SetUniform("proj_view", proj * camera->transform());

        for (auto& particle_model : particle_models)
            particle_model.draw(shader);

        for (auto& boundary_model : boundary_models)
            boundary_model.draw(shader, GL_LINES);

        shader->Unuse();

        window.SwapBuffers();
    }

public: // Helper functions
    void UpdateFrameTime()
    {
        m_frame_timer += m_elapsed_time;
        m_frame_count++;
        if (m_frame_timer >= 1.0f)
        {
            m_last_fps = m_frame_count;
            m_frame_timer -= 1.0f;
            std::println("INFO: Frame Time: {:.4f} FPS: {}", m_elapsed_time, m_frame_count);
            m_frame_count = 0;
        }
    }

};

void simulate()
{
    Simulator sim;
    if (sim.Init())
        sim.Start();
    sim.ShutDown();
}

void test()
{
    s32 n = 100;
    std::vector<f64> runtimes(n, 0.0);
    for (s32 i = 0; i < n; i++)
    {
        auto start = std::chrono::system_clock::now();

        //naive();                  // n=100 mean 41.5364ms
        //hash_map();               // n=100 mean 5.0056ms
        //unordered_map();          // n=100 mean 4.7755ms
        precompute_neighbor();   // n=100 mean 3.8802ms 

        auto end     = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration<f64, std::milli>(end - start);
        runtimes.push_back(elapsed.count());
    }
    f64 sum = 0.0;
    for (s32 i = 0; i < runtimes.size(); i++)
        sum += runtimes[i];
    f64 mean_elapsed_time = sum / n;
    std::println("n={} mean {:.4f}ms", n, mean_elapsed_time);
}

// TODO: REVIEW neighborhood search functions

#define TEST 0

int main()
{
#if TEST
    test();
#else
    simulate();
#endif
    return 0;
}