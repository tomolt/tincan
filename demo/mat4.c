#include <tincan.h>
#include <string.h>
#include <math.h>

#include "mat4.h"

#define M4(A,i,j) ((A.c)+4*(i)+(j))

Mat4
mat4_projection(tin_scalar fovy, tin_scalar aspect, tin_scalar near, tin_scalar far)
{
	Mat4 D = {{ 0 }};
	tin_scalar f = 1.0f / tanf(fovy * M_PI / 180.0f / 2.0f);
	D.c[ 0] = f / aspect;
	D.c[ 5] = f;
	D.c[10] = (near + far) / (near - far);
	D.c[11] = 2.0f * far * near / (near - far);
	D.c[14] = -1.0f;
	return D;
}

static Mat4
mat4_rotation(tin_quat q)
{
	tin_vec3 a = tin_apply_qt(q, (tin_vec3) {{ 1.0f, 0.0f, 0.0f }});
	tin_vec3 b = tin_apply_qt(q, (tin_vec3) {{ 0.0f, 1.0f, 0.0f }});
	tin_vec3 c = tin_apply_qt(q, (tin_vec3) {{ 0.0f, 0.0f, 1.0f }});
	Mat4 M = {{ 0 }};
	memcpy(M.c + 0, a.c, 3 * sizeof (float));
	memcpy(M.c + 4, b.c, 3 * sizeof (float));
	memcpy(M.c + 8, c.c, 3 * sizeof (float));
	M.c[15] = 1.0f;
	return M;
}

static Mat4
mat4_translation(tin_vec3 t)
{
	Mat4 M = {{ 0 }};
	*M4(M,3,0) = t.c[0];
	*M4(M,3,1) = t.c[1];
	*M4(M,3,2) = t.c[2];
	*M4(M,3,3) = 1.0f;
	return M;
}

static Mat4
mat4_scale(tin_scalar s)
{
	Mat4 M = {{ 0 }};
	*M4(M,0,0) = s;
	*M4(M,1,1) = s;
	*M4(M,2,2) = s;
	*M4(M,3,3) = 1.0f;
	return M;
}

Mat4
mat4_from_transform(const tin_transform *transform)
{
	Mat4 R = mat4_rotation(transform->rotation);
	Mat4 T = mat4_translation(transform->translation);
	Mat4 S = mat4_scale(transform->scale);
	return mat4_multiply(mat4_multiply(T, S), R);
}

Mat4
mat4_multiply(Mat4 A, Mat4 B)
{
	Mat4 D;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			tin_scalar f = 0.0f;
			for (int k = 0; k < 4; k++) {
				f += *M4(A,i,k) * *M4(B,k,j);
			}
			*M4(D,i,j) = f;
		}
	}
	return D;
}

