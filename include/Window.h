#pragma once

#include <print>
#include <glad/glad.h>
#include <glfw3.h>

#include "Common.h"
#include "Input.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

class Window
{
public:
     Window(std::string title, s32 width, s32 height)
	 {
         if (!glfwInit())
         {
             std::println("ERROR: Failed to initialize GLFW");
             glfwTerminate();
         }

         glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
         glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
         glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
         glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
         glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
         std::println("INFO: GLFW {}.{}.{}", GLFW_VERSION_MAJOR, GLFW_VERSION_MINOR, GLFW_VERSION_REVISION);

         window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
         if (window == nullptr)
         {
             std::println("ERROR: Failed to create GLFW window");
             glfwTerminate();
         }

         glfwMakeContextCurrent(window);

         glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
         glfwSetKeyCallback(window, key_callback);
         glfwSetMouseButtonCallback(window, mouse_button_callback);
         glfwSetCursorPosCallback(window, cursor_position_callback);
         glfwSetScrollCallback(window, scroll_callback);

         if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
         {
             std::println("ERROR: GLAD Initialization Failed");
             glfwTerminate();
         }
	 }

public:
    void PollEvents() { glfwPollEvents(); }
	bool ShouldClose() { return glfwWindowShouldClose(window); }
    void SetShouldClose() { glfwSetWindowShouldClose(window, true); }
	void Close() { glfwTerminate(); }
    void SwapBuffers() { glfwSwapBuffers(window); }
    void Clear(const ucolor uc)
    {
        fcolor fc = to_float(uc);
        glClearColor(fc.r, fc.g, fc.b, fc.a);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

private:
    GLFWwindow* window;
};