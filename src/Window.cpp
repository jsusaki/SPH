#pragma once

#include <print>
#include "Window.h"

Window::Window(std::string title, s32 width, s32 height)
{
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
    {
        std::printf("ERROR: Failed to initialize GLFW\n");
        glfwTerminate();
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_FOCUSED, GL_TRUE);
    std::printf("INFO: GLFW %d.%d.%d\n", GLFW_VERSION_MAJOR, GLFW_VERSION_MINOR, GLFW_VERSION_REVISION);

    window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    if (window == nullptr)
    {
        std::printf("ERROR: Failed to create GLFW window\n");
        glfwTerminate();
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::printf("ERROR: GLAD Initialization Failed\n");
        glfwTerminate();
    }
}

void Window::PollEvents()     { glfwPollEvents(); }
bool Window::ShouldClose()    { return glfwWindowShouldClose(window); }
void Window::SetShouldClose() { glfwSetWindowShouldClose(window, true); }
void Window::Close()          { glfwDestroyWindow(window);  glfwTerminate(); }
void Window::SwapBuffers()    { glfwSwapBuffers(window); }
void Window::Clear(const ucolor uc)
{
    fcolor fc = to_float(uc);
    glClearColor(fc.r, fc.g, fc.b, fc.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

GLFWwindow* Window::GetWindow()
{
    return window;
}
