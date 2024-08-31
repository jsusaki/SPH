// based on arcball-cpp: https://github.com/Twinklebear/arcball-cpp
#include "Common.h"

class ArcballCamera
{
public:
    ArcballCamera() = default;
    ArcballCamera(const vf3& eye, const vf3& center, const vf3& up, f32 fov, f32 aspect, f32 near, f32 far);
public:
    void init(const vf3& eye, const vf3& center, const vf3& up);

    void rotate(vf2 prev_mouse, vf2 cur_mouse);
    void pan(vf2 mouse_delta);
    void zoom(const f32 zoom_amount);
    void translate(vf3 position);

    const mf4x4& transform() const;
    const mf4x4& inv_transform() const;
    const mf4x4& projection() const;
    const mf4x4& inv_projection() const;
    const mf4x4 proj_camera() const;

    vf3 eye() const;
    vf3 dir() const;
    vf3 up() const;

private:
    void update_projection(f32 fov, f32 aspect, f32 near, f32 far);
    void update_camera();

private:
    mf4x4 center_translation;
    mf4x4 translation;
    qf32  rotation;

    mf4x4 camera;
    mf4x4 inv_camera;
    mf4x4 proj;
    mf4x4 inv_proj;

    f32 field_of_view;
    f32 aspect_ratio;
    f32 near_plane;
    f32 far_plane;
};
