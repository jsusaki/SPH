#pragma once

#include <array>

struct Input
{
    static const s32 MAX_KEYS = 257;
    static const s32 MAX_BUTTONS = 5;

    std::array<bool, MAX_KEYS>    keys      = {};
    std::array<bool, MAX_KEYS>    prev_keys = {};
    std::array<bool, MAX_BUTTONS> buttons   = {};

    vf2 mouse_pos         = {};
    vf2 prev_mouse_pos    = {};
    f32 mouse_wheel_delta = 0.0f;

    void SetKey(s32 key, bool pressed) { prev_keys[key] = keys[key]; keys[key] = pressed; }
    void SetButton(s32 button, bool pressed) { buttons[button] = pressed; }
    bool IsKeyPressed(s32 key) { return keys[key] && !prev_keys[key]; }
    bool IsKeyHeld(s32 key) { return keys[key]; }
    bool IsKeyReleased(s32 key) { return !keys[key] && prev_keys[key]; }
    bool IsButtonPressed(s32 button) { return buttons[button]; }
    void SetMousePos(f64 x, f64 y) { prev_mouse_pos = mouse_pos; mouse_pos = { static_cast<f32>(x), static_cast<f32>(y) }; }
    void SetMouseWheelDelta(f32 delta) { mouse_wheel_delta += delta; }
    void Update() { prev_keys = keys; }
};

static Input input;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS)        input.SetKey(key, true);
    else if (action == GLFW_RELEASE) input.SetKey(key, false);
}

void cursor_position_callback(GLFWwindow* window, f64 x, f64 y)
{
    input.SetMousePos(x, y);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (action == GLFW_PRESS) input.SetButton(button, true);
    else if (action == GLFW_RELEASE) input.SetButton(button, false);
}

void scroll_callback(GLFWwindow* window, f64 dx, f64 dy)
{
    for (u32 i = 0; i < std::abs(dy); i++)
    {
        if (dy > 0)      input.SetMouseWheelDelta(static_cast<f32>(1.0f));
        else if (dy < 0) input.SetMouseWheelDelta(static_cast<f32>(-1.0f));
    }
}