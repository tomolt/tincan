#include "tincan.h"

#include <tgmath.h>

/* Our matrices are column major */
#define M3(matrix,row,column) (matrix[(row)+3*(column)])

Tin_Vec3
tin_neg_v3(Tin_Vec3 x)
{
	return (Tin_Vec3) {{ -x.c[0], -x.c[1], -x.c[2] }};
}

Tin_Vec3
tin_add_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return (Tin_Vec3) {{ a.c[0]+b.c[0], a.c[1]+b.c[1], a.c[2]+b.c[2] }};
}

Tin_Vec3
tin_sub_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return (Tin_Vec3) {{ a.c[0]-b.c[0], a.c[1]-b.c[1], a.c[2]-b.c[2] }};
}

Tin_Vec3
tin_scale_v3(Tin_Scalar a, Tin_Vec3 x)
{
	return (Tin_Vec3) {{ a*x.c[0], a*x.c[1], a*x.c[2] }};
}

Tin_Vec3
tin_saxpy_v3(Tin_Scalar a, Tin_Vec3 x, Tin_Vec3 y)
{
	return tin_add_v3(tin_scale_v3(a, x), y);
}

Tin_Vec3
tin_cross_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	Tin_Vec3 c;
	c.c[0] = a.c[1] * b.c[2] - a.c[2] * b.c[1];
	c.c[1] = a.c[2] * b.c[0] - a.c[0] * b.c[2];
	c.c[2] = a.c[0] * b.c[1] - a.c[1] * b.c[0];
	return c;
}

Tin_Scalar
tin_dot_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return a.c[0]*b.c[0] + a.c[1]*b.c[1] + a.c[2]*b.c[2];
}

Tin_Scalar
tin_length_v3(Tin_Vec3 v)
{
	return sqrt(tin_dot_v3(v, v));
}

Tin_Vec3
tin_normalize_v3(Tin_Vec3 v)
{
	Tin_Scalar norm = tin_length_v3(v);
	if (norm > 0.0f) {
		return tin_scale_v3(1.0f / norm, v);
	} else {
		return v;
	}
}

Tin_Vec3
tin_hadamard_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return (Tin_Vec3) {{ a.c[0]*b.c[0], a.c[1]*b.c[1], a.c[2]*b.c[2] }};
}

void
tin_axis_angle_to_matrix(Tin_Vec3 axis, Tin_Scalar angle, Tin_Scalar matrix[3*3])
{
	Tin_Scalar cosAngle = cos(angle);
	Tin_Scalar sinAngle = sin(angle);
	Tin_Scalar antiCos = 1.0f - cosAngle;

	Tin_Scalar xyAntiCos = axis.c[0] * axis.c[1] * antiCos;
	Tin_Scalar yzAntiCos = axis.c[1] * axis.c[2] * antiCos;
	Tin_Scalar zxAntiCos = axis.c[2] * axis.c[0] * antiCos;

	Tin_Scalar xSin = axis.c[0] * sinAngle;
	Tin_Scalar ySin = axis.c[1] * sinAngle;
	Tin_Scalar zSin = axis.c[2] * sinAngle;

	matrix[0] = cosAngle + axis.c[0] * axis.c[0] * antiCos;
	matrix[1] = xyAntiCos + zSin;
	matrix[2] = zxAntiCos - ySin;
	matrix[3] = xyAntiCos - zSin;
	matrix[4] = cosAngle + axis.c[1] * axis.c[1] * antiCos;
	matrix[5] = yzAntiCos + xSin;
	matrix[6] = zxAntiCos + ySin;
	matrix[7] = yzAntiCos - xSin;
	matrix[8] = cosAngle + axis.c[2] * axis.c[2] * antiCos;
}

void
tin_m3_times_m3(Tin_Scalar result[3*3], const Tin_Scalar matrixA[3*3], const Tin_Scalar matrixB[3*3])
{
	Tin_Scalar element;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			element  = M3(matrixA,i,0) * M3(matrixB,0,j);
			element += M3(matrixA,i,1) * M3(matrixB,1,j);
			element += M3(matrixA,i,2) * M3(matrixB,2,j);
			M3(result,i,j) = element;
		}
	}
}

void
tin_m3_times_m3_transposed(Tin_Scalar result[3*3], const Tin_Scalar matrixA[3*3], const Tin_Scalar matrixB[3*3])
{
	Tin_Scalar element;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			element  = M3(matrixA,i,0) * M3(matrixB,j,0);
			element += M3(matrixA,i,1) * M3(matrixB,j,1);
			element += M3(matrixA,i,2) * M3(matrixB,j,2);
			M3(result,i,j) = element;
		}
	}
}

Tin_Vec3
tin_m3_times_v3(const Tin_Scalar matrix[3*3], Tin_Vec3 vector)
{
	Tin_Vec3 result;
	result.c[0]  = M3(matrix,0,0) * vector.c[0];
	result.c[1]  = M3(matrix,1,0) * vector.c[0];
	result.c[2]  = M3(matrix,2,0) * vector.c[0];
	result.c[0] += M3(matrix,0,1) * vector.c[1];
	result.c[1] += M3(matrix,1,1) * vector.c[1];
	result.c[2] += M3(matrix,2,1) * vector.c[1];
	result.c[0] += M3(matrix,0,2) * vector.c[2];
	result.c[1] += M3(matrix,1,2) * vector.c[2];
	result.c[2] += M3(matrix,2,2) * vector.c[2];
	return result;
}

Tin_Vec3
tin_v3_times_m3(Tin_Vec3 vector, const Tin_Scalar matrix[3*3])
{
	Tin_Vec3 result;
	for (int i = 0; i < 3; i++) {
		result.c[i]  = vector.c[0] * M3(matrix,0,i);
		result.c[i] += vector.c[1] * M3(matrix,1,i);
		result.c[i] += vector.c[2] * M3(matrix,2,i);
	}
	return result;
}

Tin_Vec3
tin_gram_schmidt(Tin_Vec3 fixed, Tin_Vec3 var)
{
	Tin_Scalar fscale = tin_dot_v3(fixed, fixed);
	if (fscale) {
		return tin_saxpy_v3(-tin_dot_v3(fixed, var) / fscale, fixed, var);
	} else {
		return var;
	}
}

Tin_Scalar
tin_prlgram_area(Tin_Vec3 e1, Tin_Vec3 e2)
{
	Tin_Vec3 perp = tin_gram_schmidt(e1, e2);
	return sqrt(tin_dot_v3(e1, e1) * tin_dot_v3(perp, perp));
}

Tin_Vec3
tin_fwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec)
{
	vec = tin_scale_v3(transform->scale, vec);
	vec = tin_fwtrf_dir(transform, vec);
	vec = tin_add_v3(transform->translation, vec);
	return vec;
}

Tin_Vec3
tin_fwtrf_dir(const Tin_Transform *transform, Tin_Vec3 vec)
{
	return tin_m3_times_v3(transform->rotation, vec);
}

Tin_Vec3
tin_bwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec)
{
	vec = tin_sub_v3(vec, transform->translation);
	vec = tin_bwtrf_dir(transform, vec);
	vec = tin_scale_v3(1.0f / transform->scale, vec);
	return vec;
}

Tin_Vec3
tin_bwtrf_dir(const Tin_Transform *transform, Tin_Vec3 vec)
{
	return tin_v3_times_m3(vec, transform->rotation);
}

