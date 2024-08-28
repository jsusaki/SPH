#include "Camera.h"

#include <cmath>
#include <iostream>

// Project the point in [-1, 1] screen space onto the arcball sphere
static qf32 screen_to_arcball(const vf2& p);

ArcballCamera::ArcballCamera(const vf3& eye, const vf3& center, const vf3& up, f32 fov, f32 aspect, f32 near, f32 far)
{
    field_of_view = fov;
    aspect_ratio  = aspect;
    near_plane    = near;
    far_plane     = far;

    init(eye, center, up);
}

void ArcballCamera::init(const vf3& eye, const vf3& center, const vf3& up)
{
    const vf3 dir = center - eye;
    vf3 z = glm::normalize(dir);
    vf3 x = glm::normalize(glm::cross(z, glm::normalize(up)));
    vf3 y = glm::normalize(glm::cross(x, z));
    x = glm::normalize(glm::cross(z, y));

    center_translation = glm::inverse(glm::translate(center));
    translation = glm::translate(vf3(0.0f, 0.0f, -glm::length(dir)));
    rotation = glm::normalize(glm::quat_cast(glm::transpose(mf3x3(x, y, -z))));

    update_camera();
}

void ArcballCamera::rotate(vf2 prev_mouse, vf2 cur_mouse)
{
    // Clamp mouse positions
    cur_mouse  = glm::clamp(cur_mouse,  vf2{ -1.0f, -1.0f }, vf2{ 1.0f, 1.0f });
    prev_mouse = glm::clamp(prev_mouse, vf2{ -1.0f, -1.0f }, vf2{ 1.0f, 1.0f });

    const qf32 mouse_cur_ball  = screen_to_arcball(cur_mouse);
    const qf32 mouse_prev_ball = screen_to_arcball(prev_mouse);

    rotation = mouse_cur_ball * mouse_prev_ball * rotation;

    update_camera();
}

void ArcballCamera::pan(vf2 mouse_delta)
{
    const f32 zoom_amount = std::abs(translation[3][2]);
    vf4 motion = { mouse_delta.x * zoom_amount, mouse_delta.y * zoom_amount, 0.0f, 0.0f };
    
    // Compute the panning in the world space
    motion = inv_camera * motion;
    center_translation = glm::translate(vf3(motion)) * center_translation;

    update_camera();
}

void ArcballCamera::zoom(const f32 zoom_amount)
{
    const vf3 motion = { 0.0f, 0.0f, zoom_amount };
    translation = glm::translate(motion) * translation;

    update_camera();
}

void ArcballCamera::translate(vf3 position)
{
    translation = glm::translate(position);

    update_camera();
}

const mf4x4& ArcballCamera::transform() const
{
    return camera;
}

const mf4x4& ArcballCamera::inv_transform() const
{
    return inv_camera;
}

vf3 ArcballCamera::eye() const
{
    return vf3{ inv_camera * vf4{0, 0, 0, 1} };
}

vf3 ArcballCamera::dir() const
{
    return glm::normalize(vf3{ inv_camera * vf4{0, 0, -1, 0} });
}

vf3 ArcballCamera::up() const
{
    return glm::normalize(vf3{ inv_camera * vf4{0, 1, 0, 0} });
}

const mf4x4& ArcballCamera::projection() const
{
    return proj;
}

const mf4x4& ArcballCamera::inv_projection() const
{
    return inv_proj;
}

const mf4x4 ArcballCamera::proj_camera() const
{
    return proj * camera;
}

void ArcballCamera::update_camera()
{
    camera      = translation * mat4_cast(rotation) * center_translation;
    inv_camera  = glm::inverse(camera);

    update_projection(field_of_view, aspect_ratio, near_plane, far_plane);
}

void ArcballCamera::update_projection(f32 fov, f32 aspect, f32 near, f32 far)
{
    proj        = glm::perspective(glm::radians(fov), aspect, near, far);
    inv_proj    = glm::inverse(proj);
}

glm::quat screen_to_arcball(const vf2& p)
{
    const float dist = glm::dot(p, p);
    // If we're on/in the sphere return the point on it
    if (dist <= 1.0f) 
    {
        return qf32(0.0f, p.x, p.y, std::sqrt(1.0f - dist));
    }
    else 
    {
        // otherwise we project the point onto the sphere
        const vf2 proj = glm::normalize(p);
        return qf32(0.0f, proj.x, proj.y, 0.f);
    }
}