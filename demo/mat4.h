typedef struct { Tin_Scalar c[16]; } Mat4;
Mat4 mat4_perspective(Tin_Scalar fovy, Tin_Scalar aspect, Tin_Scalar near, Tin_Scalar far);
Mat4 mat4_orthographic(int width, int height);
Mat4 mat4_from_transform(const Tin_Transform *transform);
Mat4 mat4_from_inverse_transform(const Tin_Transform *transform);
Mat4 mat4_multiply(Mat4 a, Mat4 b);
