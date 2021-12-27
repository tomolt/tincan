typedef struct {
	tin_quat quat;

	tin_vec3 position;
	float    yaw;
	float    pitch;
	
	bool     f, b, l, r;
	double   cursor_x;
	double   cursor_y;
} Camera;

void camera_update(Camera *camera, float inv_dt);
void camera_transform(const Camera *camera, tin_transform *transform);

