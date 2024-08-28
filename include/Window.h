#pragma once

#include <print>
#include <glad/glad.h>
#include <glfw3.h>

#include "Common.h"
#include "Input.h"
#include "Color.h"

class Window
{
public:
    Window(std::string title, s32 width, s32 height);

public:
    void PollEvents();
    bool ShouldClose();
    void SetShouldClose();
    void Close();
    void SwapBuffers();
    void Clear(const ucolor uc);

    GLFWwindow* GetWindow();

private:
    GLFWwindow* window;
};

static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

static void glfw_error_callback(int error, const char* description)
{
    std::fprintf(stderr, "ERROR: GLFW %d: %s\n", error, description);
}