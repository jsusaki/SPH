#pragma once

#include <array>
#include "Common.h"
#include <glfw3.h>

class Input
{
public: // TODO: instead of singleton make it static functions
    static Input& Instance()
    {
        static Input instance;
        return instance;
    }

    Input(const Input&) = delete;
    void operator = (const Input&) = delete;
    
    // Call this once per frame
    inline void Update() { prev_keys = keys; prev_mouse_pos = mouse_pos; }
    // Key Interface
    inline bool IsKeyPressed(const s32 key) const { return keys[key] && !prev_keys[key]; }
    inline bool IsKeyHeld(const s32 key) const { return keys[key]; }
    inline bool IsKeyReleased(const s32 key) const { return !keys[key] && prev_keys[key]; }
    // Mouse Interface
    inline vf2 GetMouse() { return mouse_pos; }
    inline bool IsButtonPressed(const s32 button) const { return buttons[button]; }
    inline f32 GetMouseWheel() { return mouse_wheel_delta; }
    inline void ResetMouseWheel() { mouse_wheel_delta = 0; }

    // Internal, used by GLFW
    inline void SetKey(const s32 key, const bool pressed) { prev_keys[key] = keys[key]; keys[key] = pressed; }
    inline void SetButton(const s32 button, const bool pressed) { buttons[button] = pressed; }
    inline void SetMousePos(const f64 x, const f64 y) { prev_mouse_pos = mouse_pos; mouse_pos = { static_cast<f32>(x), static_cast<f32>(y) }; }
    inline void SetMouseWheelDelta(const f32 delta) { mouse_wheel_delta += delta; }

private:
    Input() = default;

    static const s32 MAX_KEYS = GLFW_KEY_LAST;
    static const s32 MAX_BUTTONS = 5;
    
    std::array<bool, MAX_KEYS>    keys = {};
    std::array<bool, MAX_KEYS>    prev_keys = {};
    std::array<bool, MAX_BUTTONS> buttons = {};

    vf2 mouse_pos = {};
    vf2 prev_mouse_pos = {};
    f32 mouse_wheel_delta = 0.0f;
};

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    Input& input = Input::Instance();
    if (action == GLFW_PRESS)        input.SetKey(key, true);
    else if (action == GLFW_RELEASE) input.SetKey(key, false);
}

static void cursor_position_callback(GLFWwindow* window, f64 x, f64 y)
{
    Input& input = Input::Instance();
    input.SetMousePos(x, y);
}

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    Input& input = Input::Instance();
    if (action == GLFW_PRESS)        input.SetButton(button, true);
    else if (action == GLFW_RELEASE) input.SetButton(button, false);
}

static void scroll_callback(GLFWwindow* window, f64 dx, f64 dy)
{
    Input& input = Input::Instance();
    for (u32 i = 0; i < std::abs(dy); i++)
    {
        if (dy > 0)      input.SetMouseWheelDelta(static_cast<f32>(1.0f));
        else if (dy < 0) input.SetMouseWheelDelta(static_cast<f32>(-1.0f));
    }
}