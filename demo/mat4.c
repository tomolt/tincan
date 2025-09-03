#include <tincan.h>
#include <string.h>
#include <math.h>

#include "mat4.h"

#define M4(A,i,j) ((A.c)+4*(i)+(j))

const Mat4 mat4_identity = {{
	[ 0] = 1.0,
	[ 5] = 1.0,
	[10] = 1.0,
	[15] = 1.0,
}};

Mat4
mat4_transpose(Mat4 matrix)
{
	Mat4 transposed;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			*M4(transposed,i,j) = *M4(matrix,j,i);
		}
	}
	return transposed;
}

Mat4
mat4_perspective(Tin_Scalar fovy, Tin_Scalar aspect, Tin_Scalar near, Tin_Scalar far)
{
	Mat4 D = {{ 0 }};
	Tin_Scalar f = 1.0f / tanf(fovy * M_PI / 180.0f / 2.0f);
	*M4(D,0,0) = f / aspect;
	*M4(D,1,1) = f;
	*M4(D,2,2) = (near + far) / (near - far);
	*M4(D,2,3) = -1.0f;
	*M4(D,3,2) = 2.0f * far * near / (near - far);
	return D;
}

Mat4
mat4_orthographic(Tin_Scalar left, Tin_Scalar right, Tin_Scalar bottom, Tin_Scalar top, Tin_Scalar near, Tin_Scalar far)
{
	Mat4 D = {{ 0 }};
	*M4(D,0,0) =  2.0f / (right - left);
	*M4(D,1,1) =  2.0f / (top - bottom);
	*M4(D,2,2) = -2.0f / (far - near);
	*M4(D,3,0) = -(right + left) / (right - left);
	*M4(D,3,1) = -(top + bottom) / (top - bottom);
	*M4(D,3,2) = -(far + near) / (far - near);
	*M4(D,3,3) =  1.0f;
	return D;
}

static Mat4
mat4_rotation(const Tin_Scalar rotation[3*3])
{
	Mat4 matrix = {{ 0 }};
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix.c[4*i+j] = rotation[3*i+j];
		}
	}
	matrix.c[15] = 1.0f;
	return matrix;
}

static Mat4
mat4_translation(Tin_Vec3 t)
{
	Mat4 M = {{ 0 }};
	*M4(M,0,0) = 1.0f;
	*M4(M,1,1) = 1.0f;
	*M4(M,2,2) = 1.0f;
	*M4(M,3,0) = t.c[0];
	*M4(M,3,1) = t.c[1];
	*M4(M,3,2) = t.c[2];
	*M4(M,3,3) = 1.0f;
	return M;
}

static Mat4
mat4_scale(Tin_Scalar s)
{
	Mat4 M = {{ 0 }};
	*M4(M,0,0) = s;
	*M4(M,1,1) = s;
	*M4(M,2,2) = s;
	*M4(M,3,3) = 1.0f;
	return M;
}

Mat4
mat4_from_transform(const Tin_Transform *transform)
{
	Mat4 R = mat4_rotation(transform->rotation);
	Mat4 T = mat4_translation(transform->translation);
	Mat4 S = mat4_scale(transform->scale);
	return mat4_multiply(mat4_multiply(T, S), R);
}

Mat4
mat4_from_inverse_transform(const Tin_Transform *transform)
{
	Mat4 R = mat4_transpose(mat4_rotation(transform->rotation));
	Mat4 T = mat4_translation(tin_neg_v3(transform->translation));
	Mat4 S = mat4_scale(1.0f / transform->scale);
	return mat4_multiply(mat4_multiply(R, S), T);
}

Mat4
mat4_multiply(Mat4 A, Mat4 B)
{
	Mat4 D;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Tin_Scalar f = 0.0f;
			for (int k = 0; k < 4; k++) {
				f += *M4(A,k,j) * *M4(B,i,k);
			}
			*M4(D,i,j) = f;
		}
	}
	return D;
}

Mat4
mat4_look_at(Tin_Vec3 eye, Tin_Vec3 center, Tin_Vec3 up)
{
    Tin_Vec3 f = tin_normalize_v3(tin_sub_v3(center, eye));
    Tin_Vec3 u = tin_normalize_v3(up);
    Tin_Vec3 s = tin_normalize_v3(tin_cross_v3(f, u));
    u = tin_cross_v3(s, f);

    Mat4 M = {{ 0 }};
    *M4(M,0,0) =  s.c[0];
    *M4(M,1,0) =  s.c[1];
    *M4(M,2,0) =  s.c[2];
    *M4(M,0,1) =  u.c[0];
    *M4(M,1,1) =  u.c[1];
    *M4(M,2,1) =  u.c[2];
    *M4(M,0,2) = -f.c[0];
    *M4(M,1,2) = -f.c[1];
    *M4(M,2,2) = -f.c[2];
    *M4(M,3,0) = -tin_dot_v3(s, eye);
    *M4(M,3,1) = -tin_dot_v3(u, eye);
    *M4(M,3,2) =  tin_dot_v3(f, eye);
    *M4(M,3,3) =  1.0;
    return M;
}
