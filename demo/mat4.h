typedef struct { tin_scalar c[16]; } Mat4;
Mat4 mat4_perspective(tin_scalar fovy, tin_scalar aspect, tin_scalar near, tin_scalar far);
Mat4 mat4_orthographic(int width, int height);
Mat4 mat4_from_transform(const tin_transform *transform);
Mat4 mat4_from_inverse_transform(const tin_transform *transform);
Mat4 mat4_multiply(Mat4 a, Mat4 b);
