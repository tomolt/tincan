#include <tincan.h>
#include <stdbool.h>
#include <math.h>

#include "camera.h"

#define CAMERA_VELOCITY 3.0

void
camera_update(Camera *camera, float invDt)
{
	camera->quat = tin_make_qt((Tin_Vec3) {{ 0.0f, 1.0f, 0.0f }}, camera->yaw);
	camera->quat = tin_mul_qt(camera->quat, tin_make_qt((Tin_Vec3) {{ 1.0f, 0.0f, 0.0f }}, camera->pitch));
	
	float fb = 0.0f, lr = 0.0f;
	if (camera->f) fb += 1.0f;
	if (camera->b) fb -= 1.0f;
	if (camera->r) lr -= 1.0f;
	if (camera->l) lr += 1.0f;
	float norm = sqrtf(fb * fb + lr * lr);
	if (norm) {
		norm = (CAMERA_VELOCITY * invDt) / norm;
	}
	fb *= norm;
	lr *= norm;
	
	Tin_Vec3 fvec = tin_apply_qt(camera->quat, (Tin_Vec3) {{ 0.0f, 0.0f, -1.0f }});
	Tin_Vec3 lvec = tin_apply_qt(camera->quat, (Tin_Vec3) {{ -1.0f, 0.0f, 0.0f }});
	
	camera->position = tin_saxpy_v3(fb, fvec, camera->position);
	camera->position = tin_saxpy_v3(lr, lvec, camera->position);
}

void
camera_transform(const Camera *camera, Tin_Transform *transform)
{
	transform->translation = camera->position;
	transform->quaternion  = camera->quat;
	transform->scale       = 1.0f;
}

