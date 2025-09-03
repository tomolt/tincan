typedef struct { Tin_Scalar c[16]; } Mat4;

extern const Mat4 mat4_identity;

Mat4 mat4_perspective(Tin_Scalar fovy, Tin_Scalar aspect, Tin_Scalar near, Tin_Scalar far);
Mat4 mat4_orthographic(Tin_Scalar left, Tin_Scalar right, Tin_Scalar bottom, Tin_Scalar top, Tin_Scalar near, Tin_Scalar far);
Mat4 mat4_from_transform(const Tin_Transform *transform);
Mat4 mat4_from_inverse_transform(const Tin_Transform *transform);
Mat4 mat4_multiply(Mat4 a, Mat4 b);
Mat4 mat4_look_at(Tin_Vec3 eye, Tin_Vec3 center, Tin_Vec3 up);
