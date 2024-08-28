#include "Simulator.h"

Simulator::Simulator() : window("SPH Simulator", SCREEN_WIDTH, SCREEN_HEIGHT) 
{

}

bool Simulator::Init()
{
    // Create resources
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

bool Simulator::Start()
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

        window.SwapBuffers();

        // Update Frame Time Info
        UpdateFrameTime();
    }

    return true;
}

bool Simulator::ShutDown()
{
    gui.Shutdown();
    window.Close();
    return true;
}

void Simulator::Create()
{
    // Create Camera
    vf3 eye    = { 1.5f, 0.0f, 0.0f };
    vf3 center = { 0.0f, 0.0f, 0.0f };
    vf3 up     = { 0.0f, 1.0f, 0.0f };
    f32 fov    = 65.0f;
    f32 near   = 0.1f;
    f32 far    = 100.0f;
    f32 aspect = static_cast<f32>(SCREEN_WIDTH) / SCREEN_HEIGHT;

    camera = std::make_shared<ArcballCamera>(eye, center, up, fov, aspect, near, far);
    
    // Create SPH System Settings
    sph_settings = {
        // Simulation
        .n_particles = 25000,
        .rest_density = 997.0f,
        .gas_constant = 1.0f,
        .viscosity = 0.002f,
        .surface_tension_constant = 0.0000001f,
        .gravity = { 0.0f, -9.80665f, 0.0f },
        .mass = 0.01f,
        // Smoothing Kernel
        .smoothing_radius = 0.04f,
        .smoothing_radius2 = std::pow(0.04f, 2.0f),
        .poly6 = 315.0f / (64.0f * PI * std::pow(0.04f, 9.0f)),
        .spiky_grad = -45.0f / (PI * std::pow(0.04f, 6.0f)),
        .spiky_laplacian = 45.0f / (PI * std::pow(0.04f, 6.0f)),
        // Boundary
        .boundary_epsilon = 0.0000001f,
        .boundary_damping = -0.4f,
        .boundary_size = { 0.6f, 0.6f, 0.6f },
        .boundary_min  = { -0.6f, -0.6f, -0.6f },
        .boundary_max  = { 0.6f, 0.6f, 0.6f },
    };

    sph.Init(sph_settings);

    // Create Particle Models
    mesh sphere_mesh("./res/models/sphere.obj");
    std::vector<SPH::particle>& particles = sph.GetParticles();
    particle_models.reserve(particles.size());
    for (auto& p : particles)
    {
        model sphere_model = model(sphere_mesh);
        sphere_model.translate(p.position);
        sphere_model.scale({ 0.01f, 0.01f, 0.01f });
        particle_models.push_back(sphere_model);
    }

    // Create Bounding Box Model
    model cube_model = model(cube());
    cube_model.scale({ 
        sph_settings.boundary_size.x, 
        sph_settings.boundary_size.y, 
        sph_settings.boundary_size.z 
        });
    boundary_models.push_back(cube_model);

    // Create Shader
    shader = std::make_shared<Shader>("res/shaders/basic.vs", "res/shaders/basic.fs");

    // Create GUI
    gui.Init(window.GetWindow());
}

void Simulator::ProcessInput()
{
    Input& input = Input::Instance();

    // Key control
    if (input.IsKeyPressed(GLFW_KEY_ESCAPE))
        window.SetShouldClose();

    if (input.IsKeyPressed(GLFW_KEY_SPACE))
        simulate = !simulate;

    if (input.IsKeyHeld(GLFW_KEY_S))
        step = true;

    if (input.IsKeyPressed(GLFW_KEY_R))
        reset = true;

    if (input.IsKeyPressed(GLFW_KEY_C))
    {
        vf3 eye    = { 1.5f, 0.0f, 0.0f };
        vf3 center = { 0.0f, 0.0f, 0.0f };
        vf3 up     = { 0.0f, 1.0f, 0.0f };
        camera->init(eye, center, up);
    }

    if (input.IsKeyPressed(GLFW_KEY_TAB))
        mouse_control = !mouse_control;

    // Mouse control
    if (mouse_control)
    {
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
                vf2 dxy = mouse_pos - prev_mouse_pos;
                vf4 dxy4 = camera->inv_projection() * vf4(dxy.x, dxy.y, 0.0f, 1.0f);
                camera->pan(vf2(dxy4.x, dxy4.y));
            }
        }

        prev_mouse_pos = mouse_pos;

        if (input.mouse_wheel_delta != 0)
        {
            camera->zoom(input.mouse_wheel_delta * 0.1f);
            input.mouse_wheel_delta = 0;
        }
    }

    // Update input state
    input.Update();
    
    // Update GUI
    gui.Update();
}

void Simulator::Simulate(f32 dt)
{
    if (reset)
    {
        sph.Init(sph_settings);
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
        step = false;
    }
}

void Simulator::Render()
{
    window.Clear({ 25, 25, 25, 255 });

    // Enable rendering flags
    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    // Activate shader
    shader->Use();

    // Set Projection View Matrix
    shader->SetUniform("proj_view", camera->proj_camera());

    // Draw models
    for (auto& particle_model : particle_models)
        particle_model.draw(shader);

    for (auto& boundary_model : boundary_models)
        boundary_model.draw(shader, GL_LINES);

    // Deactivate shader
    shader->Unuse();

    // Render GUI
    gui.Render();
}

void Simulator::UpdateFrameTime()
{
    m_frame_timer += m_elapsed_time;
    m_frame_count++;
    if (m_frame_timer >= 1.0f)
    {
        m_last_fps = m_frame_count;
        m_frame_timer -= 1.0f;
        std::printf("INFO: Frame Time: %.4f FPS: %d\n", m_elapsed_time, m_frame_count);
        m_frame_count = 0;
    }
}