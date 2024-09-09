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
    camera = ArcballCamera(eye, center, up, fov, aspect, near, far);
    
    // Create SPH System Settings
    default_settings = {
        // Simulation parameters
        .n_particles              = 5000,
        .initial_particle_density = 315.0f / (64.0f * PI * std::pow(SMOOTHING_RADIUS, 3.0f)),
        .rest_density             = 1000.0f,
        .stiffness                = 3.0f,
        .viscosity                = 0.001f, //1.0f, //0.00005f,
        .surface_tension_constant = 0.0001f, //0.0728f
        .gravity                  = { 0.0f, 0.0f, 0.0f }, //{ 0.0f, -9.80665f, 0.0f },
        .mass                     = 0.02f,
        .dt                       = 0.02f,
        .radius                   = 0.03f,
        // Smoothing Kernel
        .support_radius           = 0.04f,
        .support_radius2          = std::pow(SMOOTHING_RADIUS, 2.0f),
        .poly6                    = 315.0f / (64.0f * PI * std::pow(SMOOTHING_RADIUS, 9.0f)),
        .poly6_grad               = -945.0f / (32.0f * PI * std::pow(SMOOTHING_RADIUS, 9.0f)),
        .spiky_grad               = -45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f)),
        .viscosity_laplacian      =  45.0f / (PI * std::pow(SMOOTHING_RADIUS, 6.0f)),
        // Boundary
        .boundary_epsilon = 0.0000001f,
        .boundary_damping = -1.0f,
        .boundary_size    = { 0.6f, 0.6f, 0.6f },
        .boundary_min     = {-0.6f,-0.6f,-0.6f },
        .boundary_max     = { 0.6f, 0.6f, 0.6f },
        // GFX
        .sphere_scale     = 0.0125f,
        .particle_spacing = 0.0f,
        .cube_offset      = { 0.0f, 0.25f, 0.0f },
    };
    sph.Init(default_settings);
    current_settings = default_settings;

    // Create Particle Models
    mesh sphere_mesh("./res/models/sphere.obj");
    std::vector<Particle>& particles = sph.GetParticles();
    particle_models.reserve(particles.size());
    for (auto& p : particles)
    {
        model sphere_model = model(sphere_mesh);
        sphere_model.translate(p.position);
        // gfx config
        sphere_model.scale({ 
            default_settings.sphere_scale,
            default_settings.sphere_scale,
            default_settings.sphere_scale
        });
        particle_models.push_back(sphere_model);
    }

    // Create Bounding Box Model
    model cube_model = model(cube());
    cube_model.scale({ 
        default_settings.boundary_size.x,
        default_settings.boundary_size.y,
        default_settings.boundary_size.z
    });
    boundary_models.push_back(cube_model);

    // Create Shader
    particle_shader = std::make_shared<Shader>("res/shaders/basic_particle.vert", "res/shaders/basic.frag");
    boundary_shader = std::make_shared<Shader>("res/shaders/basic_boundary.vert", "res/shaders/basic.frag");

    fcolor c = to_float(palette[0]);
    color_factor = std::vector<f32>(particles.size(), c.b);

    // Create GUI
    gui.Init(window.GetWindow());
}

void Simulator::ProcessInput()
{
    // TODO: simplify this singleton to something
    Input& input = Input::Instance();

    // Close the window
    if (input.IsKeyPressed(GLFW_KEY_ESCAPE))
        window.SetShouldClose();

    // Run the simulation
    if (input.IsKeyPressed(GLFW_KEY_SPACE))
    {
        next_state = State::SIMULATE;
        if (current_state == State::SIMULATE)
            next_state = State::PAUSE;
    }
    // Step the simulation
    if (input.IsKeyPressed(GLFW_KEY_S))
        next_state = State::STEP;
    // Restart the current simulation
    if (input.IsKeyPressed(GLFW_KEY_R) || gui.IsRestartPressed())
        next_state = State::RESTART;
    // Reset the simulation to default state
    if (input.IsKeyPressed(GLFW_KEY_T) || gui.IsResetPressed())
        next_state = State::RESET;
    // Apply the new setting to the simulation
    if (input.IsKeyPressed(GLFW_KEY_A) || gui.IsApplyPressed())
        next_state = State::APPLY;
    
    // TODO: Interpolate camera positions
    // Fixed Camera Positions
    if (input.IsKeyPressed(GLFW_KEY_1))
    {
        vf3 eye    = { 1.5f, 0.0f, 0.0f };
        vf3 center = { 0.0f, 0.0f, 0.0f };
        vf3 up     = { 0.0f, 1.0f, 0.0f };
        camera.init(eye, center, up);
    }
    if (input.IsKeyPressed(GLFW_KEY_2))
    {
        vf3 eye    = { 0.0f, 0.5f, 0.0f };
        vf3 center = { 0.0f, 0.0f, 0.0f };
        vf3 up     = { 1.0f, 0.0f, 0.0f };
        camera.init(eye, center, up);
    }
    if (input.IsKeyPressed(GLFW_KEY_3))
    {
        vf3 eye    = { 1.5f, 1.5f, 1.5f };
        vf3 center = { 0.0f, 0.0f, 0.0f };
        vf3 up     = { 0.0f, 1.0f, 0.0f };
        camera.init(eye, center, up);
    }
    if (input.IsKeyPressed(GLFW_KEY_4))
    {
        vf3 eye    = { 1.30f, 0.80f, 1.57f };
        vf3 center = { 0.0f,  0.0f, 0.0f };
        vf3 up     = { -0.22f, 0.93f,-0.29f };
        camera.init(eye, center, up);
    }
    // AVX2
    if (input.IsKeyPressed(GLFW_KEY_X)) gui.ToggleAVX2();
    sph.SetAVX2(gui.IsAVX2Enabled());
    
    // Wireframe
    if (input.IsKeyPressed(GLFW_KEY_W)) gui.ToggleWireframe();
    enable_wireframe = gui.IsWireframeEnabled();

    // Mouse control
    // TODO: incorporate mouse_pos_transformed into arcball camera
    mouse_pos_transformed = transform_mouse(input.GetMouse());
 
    if (!gui.IsWindowFocused())
    {
        if (prev_mouse_pos_transformed != (vf2(-2.0f)))
        {
            if (input.IsButtonPressed(GLFW_MOUSE_BUTTON_LEFT))
            {
                camera.rotate(prev_mouse_pos_transformed, mouse_pos_transformed);
            }
            if (input.IsButtonPressed(GLFW_MOUSE_BUTTON_RIGHT))
            {
                // TODO: incorporate into camera
                vf2 dxy  = mouse_pos_transformed - prev_mouse_pos_transformed;      // GetMouseCurrent() - GetMousePrev()
                vf4 dxy4 = camera.inv_projection() * vf4(dxy.x, dxy.y, 0.0f, 1.0f); // Get
                camera.pan(vf2(dxy4.x, dxy4.y));
            }
        }
        if (input.GetMouseWheel() != 0)
        {
            camera.zoom(input.GetMouseWheel() * 0.1f);
            input.ResetMouseWheel();
        }
    }

    prev_mouse_pos_transformed = mouse_pos_transformed;

    // Update input state
    input.Update();
}

void Simulator::Simulate(f32 dt)
{
    // Update GUI
    gui.Display(current_settings, camera, sph.GetStats());

    // Update simulation state
    switch (current_state)
    {
    case State::SIMULATE:
    case State::STEP:
    {
        sph.Simulate(current_settings);
        std::vector<Particle>& particles = sph.GetParticles();
        for (s32 i = 0; i < particles.size(); i++)
            particle_models[i].translate(particles[i].position);
        if (current_state == State::STEP)
            next_state = State::PAUSE;
        n_steps++;
        break;
    } 
    case State::RESET:
    {
        ResetSimulation(default_settings);
        current_settings = default_settings;
        next_state = State::PAUSE;
        n_steps = 0;
        break;
    } 
    case State::APPLY:
    {
        ResetSimulation(current_settings);
        next_state = State::PAUSE;
        n_steps = 0;
        break;
    } 
    case State::RESTART:
    {
        RestartSimulation();
        next_state = State::PAUSE;
        n_steps = 0;
        break;
    } 
    case State::PAUSE:
        break;
    default:
        break;
    }

    current_state = next_state;
}

void Simulator::Render()
{
    window.Clear({ 25, 25, 25, 255 });

    // Enable rendering flags
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    if (enable_wireframe)
        glEnable(GL_POLYGON_SMOOTH);
    else
        glDisable(GL_POLYGON_SMOOTH);
    
    // TODO: review  density color debugging; radius?
    /*
    if (current_state == State::SIMULATE || current_state == State::STEP)
    {
        // Get particle density
        f32 min_density = std::numeric_limits<f32>::max();
        f32 max_density = std::numeric_limits<f32>::lowest();
        // Find min and max densities
        for (s32 i = 0; i < particles.size(); i++)
        {
            Particle& p = particles[i];
            if (p.density < min_density) min_density = p.density;
            if (p.density > max_density) max_density = p.density;
        }
        // Normalize densities
        for (s32 i = 0; i < particles.size(); i++)
        {
            color_factor[i] = glm::length(p.velocity); //(p.density - min_density) / (max_density - min_density + 1e-6f);
        }
    }
    */

    // Particle Shader
    particle_shader->Use();
    particle_shader->SetUniform("proj_view", camera.proj_camera());

    f32 min_velocity = std::numeric_limits<float>::max();
    f32 max_velocity = std::numeric_limits<float>::lowest();

    std::vector<Particle>& particles = sph.GetParticles();
    for (const Particle& p : particles)
    {
        f32 velocity_magnitude = glm::length(p.velocity);
        if (velocity_magnitude < min_velocity) min_velocity = velocity_magnitude;
        if (velocity_magnitude > max_velocity) max_velocity = velocity_magnitude;
    }

    // Normalize velocity and pass as a uniform
    for (s32 i = 0; i < particle_models.size(); i++)
    {
        Particle& p = particles[i];
        f32 velocity_magnitude = glm::length(p.velocity);

        // Normalize between 0 and 1
        f32 normalized_velocity = (velocity_magnitude - min_velocity) / (max_velocity - min_velocity + 1e-6f);

        // TODO: Ability to change color palette
        // Pass normalized velocity to GPU as a color factor uniform
        particle_shader->Use();
        particle_shader->SetUniform("color_factor", normalized_velocity);
        particle_models[i].draw(particle_shader);
    }

    particle_shader->Unuse();

    // Boundary Shader
    boundary_shader->Use();
    boundary_shader->SetUniform("proj_view", camera.proj_camera());

    for (auto& boundary_model : boundary_models)
        boundary_model.draw(boundary_shader, GL_LINES);

    boundary_shader->Unuse();

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
        //std::printf("INFO: Frame Time: %.4f FPS: %d\n", m_elapsed_time, m_frame_count);
        m_frame_count = 0;
    }
}

void Simulator::ResetSimulation(const SPHSettings& settings)
{
    sph.Init(settings);
    mesh sphere_mesh("./res/models/sphere.obj");
    std::vector<Particle>& particles = sph.GetParticles();
    particle_models.clear();
    particle_models.reserve(particles.size());
    for (auto& p : particles)
    {
        model sphere_model = model(sphere_mesh);
        sphere_model.translate(p.position);
        sphere_model.scale({
            settings.sphere_scale,
            settings.sphere_scale,
            settings.sphere_scale
            });
        particle_models.push_back(sphere_model);
    }
}

void Simulator::RestartSimulation()
{
    sph.Init(current_settings);
    std::vector<Particle>& particles = sph.GetParticles();
    for (s32 i = 0; i < particles.size(); i++)
        particle_models[i].translate(particles[i].position);
}