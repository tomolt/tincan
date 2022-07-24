typedef struct {
	Tin_Quat quat;

	Tin_Vec3 position;
	float    yaw;
	float    pitch;
	
	bool     f, b, l, r;
	double   cursor_x;
	double   cursor_y;
} Camera;

void camera_update(Camera *camera, float inv_dt);
void camera_transform(const Camera *camera, Tin_Transform *transform);

