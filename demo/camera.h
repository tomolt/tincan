typedef struct {
	Tin_Scalar rotation[3*3];

	Tin_Vec3 position;
	float    yaw;
	float    pitch;
	
	bool     f, b, l, r;
	double   cursor_x;
	double   cursor_y;
} Camera;

void camera_update(Camera *camera, float invDt);
void camera_transform(const Camera *camera, Tin_Transform *transform);

