// arcball-cpp: https://github.com/Twinklebear/arcball-cpp
#include "Common.h"

class ArcballCamera
{
public:
    // Create an arcball camera focused on some center point screen: [win_width, win_height]
    ArcballCamera(const vf3& eye, const vf3& center, const vf3& up);
public:
    // Rotate the camera from the previous mouse position to the current one. Mouse positions should be in normalized device coordinates
    void rotate(vf2 prev_mouse, vf2 cur_mouse);
    // Pan the camera given the translation vector. Mouse delta amount should be in normalized device coordinates
    void pan(vf2 mouse_delta);
    // Zoom the camera given the zoom amount to (i.e., the scroll amount). Positive values zoom in, negative will zoom out.
    void zoom(const f32 zoom_amount);
    // Get the camera transformation matrix
    const mf4x4& transform() const;
    // Get the camera's inverse transformation matrix
    const mf4x4& inv_transform() const;
    // Get the eye position of the camera in world space
    vf3 eye() const;
    // Get the eye direction of the camera in world space
    vf3 dir() const;
    // Get the up direction of the camera in world space
    vf3 up() const;
private:
    void update_camera();
private:
    // We store the unmodified look at matrix along with decomposed translation and rotation components
    mf4x4 center_translation;
    mf4x4 translation;
    qf32  rotation;
    // Camera is the full camera transform, inv_camera is stored as well to easily compute eye position and world space rotation axes
    mf4x4 camera;
    mf4x4 inv_camera;

    // projection matrix
    // view matrix
};
