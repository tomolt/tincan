#include <tincan.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "camera.h"

#define CAMERA_VELOCITY 3.0

void
camera_update(Camera *camera, float invDt)
{
	Tin_Scalar yawMatrix[3*3];
	tin_axis_angle_to_matrix((Tin_Vec3){{ 0.0f, 1.0f, 0.0f }}, camera->yaw, yawMatrix);

	Tin_Scalar pitchMatrix[3*3];
	tin_axis_angle_to_matrix((Tin_Vec3){{ 1.0f, 0.0f, 0.0f }}, camera->pitch, pitchMatrix);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			camera->rotation[i+3*j]  = yawMatrix[i+3*0] * pitchMatrix[0+3*j];
			camera->rotation[i+3*j] += yawMatrix[i+3*1] * pitchMatrix[1+3*j];
			camera->rotation[i+3*j] += yawMatrix[i+3*2] * pitchMatrix[2+3*j];
		}
	}

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

	Tin_Vec3 fvec = {{ -camera->rotation[6], -camera->rotation[7], -camera->rotation[8] }};
	Tin_Vec3 lvec = {{ -camera->rotation[0], -camera->rotation[1], -camera->rotation[2] }};
	
	camera->position = tin_saxpy_v3(fb, fvec, camera->position);
	camera->position = tin_saxpy_v3(lr, lvec, camera->position);
}

void
camera_transform(const Camera *camera, Tin_Transform *transform)
{
	memcpy(transform->rotation, camera->rotation, sizeof camera->rotation);
	transform->translation = camera->position;
	transform->scale = 1.0f;
}

