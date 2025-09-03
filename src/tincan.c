/* tincan physics engine
 *
 * Copyright 2021-2025 Thomas Oltmann
 * All rights reserved
 */

#include "tincan.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <assert.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* === Vector Math === :vec: */

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
	for (int i = 0; i < 3; i++) {
		result.c[i]  = M3(matrix,i,0) * vector.c[0];
		result.c[i] += M3(matrix,i,1) * vector.c[1];
		result.c[i] += M3(matrix,i,2) * vector.c[2];
	}
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

/* === Transforms === :trf: */

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

/* === Polytopes === :poly: */

Tin_Vec3
tin_polytope_support(const void *geometry, Tin_Vec3 dir)
{
	const Tin_Polytope *polytope = geometry;
	Tin_Scalar bestScore = -INFINITY;
	int bestIdx = -1;
	for (int idx = 0; idx < polytope->numVertices; idx++) {
		Tin_Scalar score = tin_dot_v3(polytope->vertices[idx], dir);
		if (score > bestScore) {
			bestScore = score;
			bestIdx   = idx;
		}
	}
	return polytope->vertices[bestIdx];
}

Tin_Vec3
tin_polysum_support(const void *geometry, Tin_Vec3 dir)
{
	const Tin_Polysum *s = geometry;

	Tin_Vec3 former_dir = tin_bwtrf_dir(s->transform1, dir);
	Tin_Vec3 relTo1 = tin_polytope_support(s->polytope1, former_dir);
	Tin_Vec3 former_abs = tin_fwtrf_point(s->transform1, relTo1);

	Tin_Vec3 latter_dir = tin_bwtrf_dir(s->transform2, tin_neg_v3(dir));
	Tin_Vec3 relTo2 = tin_polytope_support(s->polytope2, latter_dir);
	Tin_Vec3 latter_abs = tin_fwtrf_point(s->transform2, relTo2);

	return tin_sub_v3(former_abs, latter_abs);
}

/* === Shapes === :shape: */

void
tin_shape_aabb(const Tin_Shape *shape, const Tin_Transform *transform, Tin_Vec3 *aabbMin, Tin_Vec3 *aabbMax)
{
	switch (shape->kind) {
	case TIN_SPHERE:
		*aabbMin = (Tin_Vec3){{ -1.0, -1.0, -1.0 }};
		*aabbMax = (Tin_Vec3){{  1.0,  1.0,  1.0 }};
		*aabbMin = tin_saxpy_v3(transform->scale, *aabbMin, transform->translation);
		*aabbMax = tin_saxpy_v3(transform->scale, *aabbMax, transform->translation);
		break;

	case TIN_POLYTOPE:
		// TODO reuse back-transformed cardinal direction vectors
		aabbMin->c[0] = tin_polytope_support(&shape->polytope, tin_bwtrf_dir(transform, (Tin_Vec3){{ -1.0, 0.0, 0.0 }})).c[0];
		aabbMin->c[1] = tin_polytope_support(&shape->polytope, tin_bwtrf_dir(transform, (Tin_Vec3){{  0.0,-1.0, 0.0 }})).c[1];
		aabbMin->c[2] = tin_polytope_support(&shape->polytope, tin_bwtrf_dir(transform, (Tin_Vec3){{  0.0, 0.0,-1.0 }})).c[2];
		aabbMax->c[0] = tin_polytope_support(&shape->polytope, tin_bwtrf_dir(transform, (Tin_Vec3){{  1.0, 0.0, 0.0 }})).c[0];
		aabbMax->c[1] = tin_polytope_support(&shape->polytope, tin_bwtrf_dir(transform, (Tin_Vec3){{  0.0, 1.0, 0.0 }})).c[1];
		aabbMax->c[2] = tin_polytope_support(&shape->polytope, tin_bwtrf_dir(transform, (Tin_Vec3){{  0.0, 0.0, 1.0 }})).c[2];
		*aabbMin = tin_fwtrf_point(transform, *aabbMin);
		*aabbMax = tin_fwtrf_point(transform, *aabbMax);
		break;
	
	default:
		abort();
	}
}

/* === Minkowski Portal Refinement === :mpr: */

int
tin_construct_portal(const void *geometry, Tin_SupportFunc support,
	const Tin_Ray *r, Tin_Portal *p)
{
	/* find the very first point */
	p->a = support(geometry, r->dir);

	/* find the second point */
	{
		Tin_Vec3 oa = tin_sub_v3(p->a, r->origin);
		Tin_Vec3 dir = tin_gram_schmidt(oa, r->dir);
		if (tin_dot_v3(dir, dir) == 0.0f) {
			dir = tin_gram_schmidt(oa, (Tin_Vec3) {{ 1.0f, 0.0f, 0.0f }});
			if (tin_dot_v3(dir, dir) == 0.0f) {
				dir = (Tin_Vec3) {{ 0.0f, 0.0f, 1.0f }};
			}
		}
		dir = tin_normalize_v3(dir);
		p->b = support(geometry, dir);
	}

	for (int it = 0;; it++) {
		if (it >= 100) {
			fprintf(stderr, "Portal construction took too many iterations.\n");
			return 0;
		}

		/*  */
		Tin_Vec3 oa = tin_sub_v3(p->a, r->origin);
		Tin_Vec3 ob = tin_sub_v3(p->b, r->origin);
		Tin_Vec3 dir = tin_cross_v3(oa, ob);
		if (tin_dot_v3(dir, r->dir) < 0.0f) {
			dir = tin_neg_v3(dir);
		}
		p->c = support(geometry, dir);
		Tin_Scalar score = tin_dot_v3(dir, p->c);
		if (score <= tin_dot_v3(dir, p->a) || score <= tin_dot_v3(dir, p->b)) {
			return 0;
		}
		Tin_Vec3 oc = tin_sub_v3(p->c, r->origin);
		
		/*  */
		Tin_Vec3 an = tin_cross_v3(ob, oc);
		if (tin_dot_v3(an, oa) * tin_dot_v3(an, r->dir) < 0.0f) {
			p->a = p->c;
			continue;
		}
		
		/*  */
		Tin_Vec3 bn = tin_cross_v3(oc, oa);
		if (tin_dot_v3(bn, ob) * tin_dot_v3(bn, r->dir) < 0.0f) {
			p->b = p->c;
			continue;
		}
		
		return 1;
	}
}

void
tin_refine_portal(const void *geometry, Tin_SupportFunc support,
	const Tin_Ray *r, Tin_Portal *p)
{
	Tin_Vec3 ab = tin_sub_v3(p->b, p->a);
	Tin_Vec3 ac = tin_sub_v3(p->c, p->a);
	p->normal = tin_cross_v3(ab, ac);
	if (tin_dot_v3(p->normal, r->dir) < 0.0f) {
		Tin_Vec3 temp;
		temp = p->a;
		p->a = p->b;
		p->b = temp;
	}

	for (int it = 0;; it++) {
		Tin_Vec3 ab = tin_sub_v3(p->b, p->a);
		Tin_Vec3 ac = tin_sub_v3(p->c, p->a);
		p->normal = tin_cross_v3(ab, ac);

		if (tin_dot_v3(p->normal, p->normal) == 0.0f) {
			fprintf(stderr, "portal collapsed. (it=%d)\n", it);
			break;
		}

		Tin_Vec3 s = support(geometry, p->normal);
		
		{
			Tin_Scalar sScore = tin_dot_v3(s, p->normal);
			if (sScore <= 0.0f) break;
			Tin_Scalar aScore = tin_dot_v3(p->a, p->normal);
			if (sScore <= aScore * 1.01f) break;
			if (sScore - aScore <= 1e-6f) break;
		}
		
		if (it >= 100) {
			fprintf(stderr, "MPR refine_portal() took too many iterations.\n");
			fprintf(stderr, "\tnormal length squared: %f\n", tin_dot_v3(p->normal, p->normal));
			Tin_Scalar aScore = tin_dot_v3(p->a, p->normal);
			Tin_Scalar sScore = tin_dot_v3(s, p->normal);
			fprintf(stderr, "\taScore: %e\n", aScore);
			fprintf(stderr, "\tsScore: %e\n", sScore);
			break;
		}

		Tin_Vec3 os = tin_sub_v3(s, r->origin);
		Tin_Vec3 d_x_os = tin_cross_v3(r->dir, os);
		
		/* if <D , (OS x OA)> > 0 */
		Tin_Vec3 oa = tin_sub_v3(p->a, r->origin);
		if (tin_dot_v3(oa, d_x_os) > 0.0f) {
			/* if <D , (OS x OB)> > 0 */
			Tin_Vec3 ob = tin_sub_v3(p->b, r->origin);
			if (tin_dot_v3(ob, d_x_os) > 0.0f) {
				p->a = s; /* SBC */
			} else {
				p->c = s; /* ABS */
			}
		} else {
			/* if <D , (OS x OC)> > 0 */
			Tin_Vec3 oc = tin_sub_v3(p->c, r->origin);
			if (tin_dot_v3(oc, d_x_os) > 0.0f) {
				p->b = s; /* ASC */
			} else {
				p->a = s; /* SBC */
			}
		}
	}
}

/* === Contact Points === :contact: */

int
tin_clip_manifold_against_plane(const Tin_Vec3 *points, int count, Tin_Vec3 normal, Tin_Scalar base, Tin_Vec3 *newPoints)
{
	int newCount = 0;
	int j = count - 1;
	Tin_Scalar pj = tin_dot_v3(normal, points[j]) - base;
	for (int i = 0; i < count; i++) {
		Tin_Scalar pi = tin_dot_v3(normal, points[i]) - base;
		if (pi <= 0.0) {
			if (!(pj <= 0.0)) {
				Tin_Vec3 dir = tin_sub_v3(points[i], points[j]);
				Tin_Scalar t = (base - tin_dot_v3(normal, points[j])) / tin_dot_v3(normal, dir);
				Tin_Vec3 x = tin_saxpy_v3(t, dir, points[j]);
				newPoints[newCount++] = x;
			}
			newPoints[newCount++] = points[i];
		} else {
			if (pj <= 0.0) {
				Tin_Vec3 dir = tin_sub_v3(points[i], points[j]);
				Tin_Scalar t = (base - tin_dot_v3(normal, points[j])) / tin_dot_v3(normal, dir);
				Tin_Vec3 x = tin_saxpy_v3(t, dir, points[j]);
				newPoints[newCount++] = x;
			}
		}
		j = i;
		pj = pi;
	}
	return newCount;
}

int
tin_clip_manifolds(const Tin_Vec3 *pointsA, int countA, const Tin_Vec3 *pointsB, int countB, Tin_Vec3 perpendicular, Tin_Vec3 *newPoints)
{
	int newCount = countB;
	memcpy(newPoints, pointsB, newCount * sizeof *newPoints);
	int j = countA - 1;
	for (int i = 0; i < countA; i++) {
		if (newCount == 0) break;
		Tin_Vec3 edge = tin_sub_v3(pointsA[i], pointsA[j]);
		Tin_Vec3 normal = tin_cross_v3(edge, perpendicular); // TODO determine the correct winding here
		normal = tin_normalize_v3(normal);
		Tin_Scalar base = tin_dot_v3(normal, pointsA[i]);
		Tin_Vec3 buffer[32];
		printf("%d -> ", newCount);
		newCount = tin_clip_manifold_against_plane(newPoints, newCount, normal, base, buffer);
		printf("%d\n", newCount);
		memcpy(newPoints, buffer, newCount * sizeof *newPoints);
		j = i;
	}
	return newCount;
}

int
tin_incident_face(const Tin_Polytope *polytope, Tin_Vec3 dir)
{
	int bestFace = -1;
	Tin_Scalar bestScore = -INFINITY;
	for (int f = 0; f < polytope->numFaces; f++) {
		Tin_Scalar score = tin_dot_v3(dir, polytope->faceNormals[f]);
		if (score > bestScore) {
			bestFace = f;
			bestScore = score;
		}
	}
	return bestFace;
}

int
tin_reduce_manifold(Tin_Vec3 *points, int count)
{
	assert(count > 0);
	int idx = -1;
	Tin_Scalar bestScore = INFINITY;
	for (int i = 0; i < count; i++) {
		int j1 = i - 1;
		if (j1 < 0) j1 = count - 1;
		int j2 = i + 1;
		if (j2 >= count) j2 = 0;
		Tin_Vec3 e1 = tin_sub_v3(points[i], points[j1]);
		Tin_Vec3 e2 = tin_sub_v3(points[j2], points[i]);
		Tin_Scalar score = tin_prlgram_area(e1, e2);
		if (score < bestScore) {
			idx = i;
			bestScore = score;
		}
	}
	assert(idx >= 0);
	assert(idx < count);
	count -= 1;
	if (count - idx > 0) {
		memmove(points + idx, points + idx + 1, (count - idx) * sizeof *points);
	}
	return count;
}

int
tin_polytope_collide(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Contact *contacts)
{
	Tin_Polysum ps = { pa, ta, pb, tb };

	Tin_Ray r;
	r.origin = tin_sub_v3(ta->translation, tb->translation);
	r.dir    = tin_neg_v3(r.origin);
	Tin_Scalar norm = sqrt(tin_dot_v3(r.dir, r.dir));
	if (norm == 0.0f) {
		/* FIXME */
		return 0;
	}
	r.dir = tin_normalize_v3(r.dir);

	Tin_Portal p;
	if (!tin_construct_portal(&ps, tin_polysum_support, &r, &p)) {
		return 0;
	}
	tin_refine_portal(&ps, tin_polysum_support, &r, &p);
	p.normal = tin_normalize_v3(p.normal);
	if (tin_dot_v3(p.normal, p.a) <= 0.0f) {
		return 0;
	}

	for (int it = 0; it < 4; it++) {
		Tin_Ray nr;
		nr.dir    = p.normal;
		nr.origin = (Tin_Vec3) {{ 0.0f, 0.0f, 0.0f }};
		Tin_Portal np;
		if (!tin_construct_portal(&ps, tin_polysum_support, &nr, &np)) {
			break;
		}
		tin_refine_portal(&ps, tin_polysum_support, &nr, &np);
		np.normal = tin_normalize_v3(np.normal);

		Tin_Scalar proj = tin_dot_v3(p.normal, np.normal);
		p = np; /* lol */
		r = nr;
		if (proj >= 0.99f) break;
	}

	int faceA = tin_incident_face(pa, tin_bwtrf_dir(ta, p.normal));
	printf("faceA: %d\n", faceA);
	Tin_Vec3 pointsA[32];
	int countA = 0;
	for (int i = pa->faceOffsets[faceA]; i < pa->faceOffsets[faceA+1]; i++) {
		int idx = pa->faceIndices[i];
		pointsA[countA++] = tin_fwtrf_point(ta, pa->vertices[idx]);
	}

	Tin_Vec3 refNormal = tin_fwtrf_dir(ta, pa->faceNormals[faceA]);
	Tin_Scalar refBase = tin_dot_v3(refNormal, pointsA[0]);

	int faceB = tin_incident_face(pb, tin_bwtrf_dir(tb, tin_neg_v3(refNormal)));
	Tin_Vec3 pointsB[32];
	int countB = 0;
	for (int i = pb->faceOffsets[faceB]; i < pb->faceOffsets[faceB+1]; i++) {
		int idx = pb->faceIndices[i];
		pointsB[countB++] = tin_fwtrf_point(tb, pb->vertices[idx]);
	}

	Tin_Vec3 manifold[32];
	int count = tin_clip_manifolds(pointsA, countA, pointsB, countB, refNormal, manifold);

	while (count > TIN_MAX_CONTACTS) {
		count = tin_reduce_manifold(manifold, count);
	}

	memset(contacts, 0, count * sizeof *contacts);
	for (int i = 0; i < count; i++) {
		Tin_Scalar projDist = tin_dot_v3(refNormal, manifold[i]) - refBase;
		Tin_Vec3 projPoint = tin_saxpy_v3(-projDist, refNormal, manifold[i]);
		contacts[i].rel1 = tin_bwtrf_point(ta, projPoint);
		contacts[i].rel2 = tin_bwtrf_point(tb, manifold[i]);
		contacts[i].normal = refNormal;
		contacts[i].rel1Normal = tin_bwtrf_dir(ta, refNormal);
		contacts[i].position = manifold[i];
	}

	return count;
}

Tin_Scalar
tin_dot_array(Tin_Scalar a[], Tin_Scalar b[], int length)
{
	Tin_Scalar sum = 0.0f;
	for (int i = 0; i < length; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

Tin_Vec3
tin_solve_inertia(const Tin_Body *body, Tin_Vec3 vec)
{
	return tin_m3_times_v3(body->invInertia, vec);
}

Tin_Scalar
tin_effective_mass(const Tin_Body *body1, const Tin_Body *body2, Tin_Scalar jacobian[12])
{
	Tin_Scalar invMassJacobian[12];
	Tin_Vec3 temp;
	invMassJacobian[0] = body1->invMass * jacobian[0];
	invMassJacobian[1] = body1->invMass * jacobian[1];
	invMassJacobian[2] = body1->invMass * jacobian[2];
	temp = tin_solve_inertia(body1, (Tin_Vec3){{ jacobian[3], jacobian[4], jacobian[5] }});
	invMassJacobian[3] = temp.c[0];
	invMassJacobian[4] = temp.c[1];
	invMassJacobian[5] = temp.c[2];
	invMassJacobian[6] = body2->invMass * jacobian[6];
	invMassJacobian[7] = body2->invMass * jacobian[7];
	invMassJacobian[8] = body2->invMass * jacobian[8];
	temp = tin_solve_inertia(body2, (Tin_Vec3){{ jacobian[9], jacobian[10], jacobian[11] }});
	invMassJacobian[9] = temp.c[0];
	invMassJacobian[10] = temp.c[1];
	invMassJacobian[11] = temp.c[2];

	Tin_Scalar effectiveMass = 1.0f / tin_dot_array(jacobian, invMassJacobian, 12);
	return effectiveMass;
}

void
tin_apply_impulse(Tin_Body *body1, Tin_Body *body2, Tin_Scalar jacobian[12], Tin_Scalar magnitude)
{
	body1->velocity = tin_saxpy_v3(magnitude * body1->invMass, TIN_VEC3(jacobian[0], jacobian[1], jacobian[2]), body1->velocity);
	body1->angularVelocity = tin_saxpy_v3(magnitude, tin_solve_inertia(body1, TIN_VEC3(jacobian[3], jacobian[4], jacobian[5])), body1->angularVelocity);
	body2->velocity = tin_saxpy_v3(magnitude * body2->invMass, TIN_VEC3(jacobian[6], jacobian[7], jacobian[8]), body2->velocity);
	body2->angularVelocity = tin_saxpy_v3(magnitude, tin_solve_inertia(body2, TIN_VEC3(jacobian[9], jacobian[10], jacobian[11])), body2->angularVelocity);
}

void
tin_enforce_jacobian(Tin_Body *body1, Tin_Body *body2, Tin_Scalar jacobian[12],
	Tin_Scalar effectiveMass, Tin_Scalar bias,
	Tin_Scalar *magnitudeAccum, Tin_Scalar minMagnitude, Tin_Scalar maxMagnitude)
{
	/* Solve for magnitude */
	Tin_Scalar magnitude = bias;
	magnitude -= tin_dot_array(jacobian + 0, body1->velocity.c, 3);
	magnitude -= tin_dot_array(jacobian + 3, body1->angularVelocity.c, 3);
	magnitude -= tin_dot_array(jacobian + 6, body2->velocity.c, 3);
	magnitude -= tin_dot_array(jacobian + 9, body2->angularVelocity.c, 3);
	magnitude *= effectiveMass;

	/* Clamp magnitude to fulfill inequality */
	if (magnitudeAccum != NULL) {
		Tin_Scalar prevAccum = *magnitudeAccum;
		*magnitudeAccum += magnitude;
		*magnitudeAccum = MAX(*magnitudeAccum, minMagnitude);
		*magnitudeAccum = MIN(*magnitudeAccum, maxMagnitude);
		magnitude = *magnitudeAccum - prevAccum;
	}

	tin_apply_impulse(body1, body2, jacobian, magnitude);
}

void
tin_arbiter_add_contact(Tin_Arbiter *arbiter, Tin_Contact contact)
{
	Tin_Vec3 p1, p2;
	p1 = tin_fwtrf_point(&arbiter->body1->transform, contact.rel1);
	p2 = tin_fwtrf_point(&arbiter->body2->transform, contact.rel2);
	contact.baseStretch = tin_length_v3(tin_gram_schmidt(contact.normal, tin_sub_v3(p1, p2)));

	Tin_Vec3 newPos = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
	int minIdx = TIN_MAX_CONTACTS;
	Tin_Scalar minDist = INFINITY;
	for (int idx = 0; idx < arbiter->numContacts; idx++) {
		Tin_Contact *c = &arbiter->contacts[idx];
		p1 = tin_fwtrf_point(&arbiter->body1->transform, c->rel1);
		p2 = tin_fwtrf_point(&arbiter->body2->transform, c->rel2);
		Tin_Vec3 oldPos = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
		Tin_Vec3 diff = tin_sub_v3(oldPos, newPos);
		Tin_Scalar dist = tin_dot_v3(diff, diff);
		if (dist < minDist) {
			minIdx  = idx;
			minDist = dist;
		}
	}

	/*if (minDist < 0.001) {
		return;
	}*/

	int newIdx;
	if (arbiter->numContacts < TIN_MAX_CONTACTS) {
		newIdx = arbiter->numContacts++;
	} else {
		newIdx = minIdx;
	}

	contact.ineqAccum = 0.0f;
	contact.tangentAccum = 0.0f;
	contact.bitangentAccum = 0.0f;
	arbiter->contacts[newIdx] = contact;
}

void
tin_jacobian_along_axis(Tin_Scalar jacobian[12], Tin_Vec3 axis, Tin_Vec3 r1, Tin_Vec3 r2)
{
	Tin_Vec3 r1Xaxis = tin_cross_v3(r1, axis);
	Tin_Vec3 r2Xaxis = tin_cross_v3(r2, axis);
	jacobian[0] = -axis.c[0];
	jacobian[1] = -axis.c[1];
	jacobian[2] = -axis.c[2];
	jacobian[3] = -r1Xaxis.c[0];
	jacobian[4] = -r1Xaxis.c[1];
	jacobian[5] = -r1Xaxis.c[2];
	jacobian[6] = axis.c[0];
	jacobian[7] = axis.c[1];
	jacobian[8] = axis.c[2];
	jacobian[9] = r2Xaxis.c[0];
	jacobian[10] = r2Xaxis.c[1];
	jacobian[11] = r2Xaxis.c[2];
}

void
tin_arbiter_update(Tin_Arbiter *arbiter)
{
	const Tin_Scalar maxSeparation = 0.1f;
	const Tin_Scalar maxStretch = 0.2f;

	for (int i = 0; i < arbiter->numContacts; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		Tin_Vec3 p1 = tin_fwtrf_point(&arbiter->body1->transform, contact->rel1);
		Tin_Vec3 p2 = tin_fwtrf_point(&arbiter->body2->transform, contact->rel2);

		Tin_Scalar separation = tin_dot_v3(contact->normal, tin_sub_v3(p2, p1));

		if (separation > maxSeparation) {
			*contact = arbiter->contacts[--arbiter->numContacts];
			i--;
			continue;
		}

		Tin_Scalar stretch = tin_length_v3(tin_gram_schmidt(contact->normal, tin_sub_v3(p1, p2)));
		if (stretch > maxStretch) {
			*contact = arbiter->contacts[--arbiter->numContacts];
			i--;
			continue;
		}
	}

}

void
tin_arbiter_prestep(Tin_Arbiter *arbiter, Tin_Scalar invDt)
{
	const Tin_Scalar allowedPenetration = 0.005;
	const Tin_Scalar biasFactor = 0.1;

	for (int i = 0; i < arbiter->numContacts; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		Tin_Vec3 p1 = tin_fwtrf_point(&arbiter->body1->transform, contact->rel1);
		Tin_Vec3 p2 = tin_fwtrf_point(&arbiter->body2->transform, contact->rel2);

		//contact->position   = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
		contact->position = p1;
		contact->normal = tin_fwtrf_dir(&arbiter->body1->transform, contact->rel1Normal);
		contact->separation = tin_dot_v3(contact->normal, tin_sub_v3(p2, p1));

		Tin_Vec3 r1 = tin_sub_v3(contact->position, arbiter->body1->transform.translation);
		Tin_Vec3 r2 = tin_sub_v3(contact->position, arbiter->body2->transform.translation);

		/* Precompute jacobian, effectiveMass, and bias */
		tin_jacobian_along_axis(contact->jacobian, contact->normal, r1, r2);
		contact->effectiveMass[0] = tin_effective_mass(arbiter->body1, arbiter->body2, contact->jacobian);

		contact->bias = -biasFactor * invDt * MIN(0.0f, contact->separation + allowedPenetration) / arbiter->numContacts;

		contact->ineqAccum = 0.0f;
		contact->tangentAccum = 0.0f;
		contact->bitangentAccum = 0.0f;
	}
}

void
tin_arbiter_warm_start(Tin_Arbiter *arbiter)
{
	for (int idx = 0; idx < arbiter->numContacts; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];
		if (contact->separation >= 0.0f) continue;

		tin_apply_impulse(arbiter->body1, arbiter->body2, contact->jacobian, contact->ineqAccum);
	}
}

void
tin_arbiter_apply_separation(Tin_Arbiter *arbiter, Tin_Scalar invDt)
{
	// TODO
	(void)invDt;

	for (int idx = 0; idx < arbiter->numContacts; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];
		if (contact->separation >= 0.0f) continue;

		tin_enforce_jacobian(arbiter->body1, arbiter->body2, contact->jacobian, contact->effectiveMass[0], contact->bias, &contact->ineqAccum, 0.0f, INFINITY);
	}
}

void
tin_arbiter_apply_friction(Tin_Arbiter *arbiter, Tin_Scalar invDt)
{
	// TODO
	(void)invDt;

	const Tin_Scalar friction = 0.5f;

	for (int idx = 0; idx < arbiter->numContacts; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];
		if (contact->separation >= 0.0f) continue;

		Tin_Vec3 r1 = tin_fwtrf_dir(&arbiter->body1->transform,
			tin_scale_v3(arbiter->body1->transform.scale, contact->rel1));
		Tin_Vec3 r2 = tin_fwtrf_dir(&arbiter->body2->transform,
			tin_scale_v3(arbiter->body2->transform.scale, contact->rel2));

		Tin_Vec3 tangent1 = tin_gram_schmidt(contact->normal, (Tin_Vec3){{ 1.0f, 0.0f, 0.0f }});
		if (tin_dot_v3(tangent1, tangent1) == 0.0f) {
			tangent1 = tin_gram_schmidt(contact->normal, (Tin_Vec3){{ 0.0f, 0.0f, 1.0f }});
			if (tin_dot_v3(tangent1, tangent1) == 0.0f) {
				tangent1 = (Tin_Vec3){{ 0.0f, 1.0f, 0.0f }};
			}
		}
		Tin_Vec3 tangent2 = tin_cross_v3(contact->normal, tangent1);

		Tin_Scalar jacobian[12];
		Tin_Scalar effectiveMass;

		tin_jacobian_along_axis(jacobian, tangent1, r1, r2);
		effectiveMass = tin_effective_mass(arbiter->body1, arbiter->body2, jacobian);
		tin_enforce_jacobian(arbiter->body1, arbiter->body2, jacobian, effectiveMass, 0.0f, &contact->tangentAccum, -contact->ineqAccum * friction, contact->ineqAccum * friction);

		tin_jacobian_along_axis(jacobian, tangent2, r1, r2);
		effectiveMass = tin_effective_mass(arbiter->body1, arbiter->body2, jacobian);
		tin_enforce_jacobian(arbiter->body1, arbiter->body2, jacobian, effectiveMass, 0.0f, &contact->bitangentAccum, -contact->ineqAccum * friction, contact->ineqAccum * friction);
	}
}

void
tin_joint_apply_impulse(Tin_Joint *joint, Tin_Scalar invDt)
{
	const Tin_Scalar biasFactor = 0.1f;

	Tin_Scalar bias;
	Tin_Scalar separation;
	Tin_Scalar jacobian[12];
	Tin_Scalar effectiveMass;

	Tin_Vec3 p1 = tin_fwtrf_point(&joint->body1->transform, joint->relTo1);
	Tin_Vec3 p2 = tin_fwtrf_point(&joint->body2->transform, joint->relTo2);

	Tin_Vec3 r1 = tin_fwtrf_dir(&joint->body1->transform,
		tin_scale_v3(joint->body1->transform.scale, joint->relTo1));
	Tin_Vec3 r2 = tin_fwtrf_dir(&joint->body2->transform,
		tin_scale_v3(joint->body2->transform.scale, joint->relTo2));

	separation = p2.c[0] - p1.c[0];
	bias = -biasFactor * invDt * separation;
	tin_jacobian_along_axis(jacobian, (Tin_Vec3){{1.0f, 0.0f, 0.0f}}, r1, r2);
	effectiveMass = tin_effective_mass(joint->body1, joint->body2, jacobian);
	tin_enforce_jacobian(joint->body1, joint->body2, jacobian, effectiveMass, bias, NULL, 0.0f, 0.0f);

	separation = p2.c[1] - p1.c[1];
	bias = -biasFactor * invDt * separation;
	tin_jacobian_along_axis(jacobian, (Tin_Vec3){{0.0f, 1.0f, 0.0f}}, r1, r2);
	effectiveMass = tin_effective_mass(joint->body1, joint->body2, jacobian);
	tin_enforce_jacobian(joint->body1, joint->body2, jacobian, effectiveMass, bias, NULL, 0.0f, 0.0f);

	separation = p2.c[2] - p1.c[2];
	bias = -biasFactor * invDt * separation;
	tin_jacobian_along_axis(jacobian, (Tin_Vec3){{0.0f, 0.0f, 1.0f}}, r1, r2);
	effectiveMass = tin_effective_mass(joint->body1, joint->body2, jacobian);
	tin_enforce_jacobian(joint->body1, joint->body2, jacobian, effectiveMass, bias, NULL, 0.0f, 0.0f);
}

/* === Pair-Indexed Hashtable === :pair: */

void
tin_create_pairtable(Tin_PairTable *table)
{
	table->count = 0;
	table->capac = 16;
	table->slots = calloc(table->capac, sizeof *table->slots);
}

void
tin_destroy_pairtable(Tin_PairTable *table)
{
	free(table->slots);
}

void
tin_reset_pairtable(Tin_PairTable *table)
{
	table->count = 0;
	memset(table->slots, 0, table->capac * sizeof *table->slots);
}

/* Capacity has to be a power of two.
 */
void
tin_resize_pairtable(Tin_PairTable *table, size_t newCapac)
{
	void tin_insert_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void *payload);

	size_t oldCapac = table->capac;
	Tin_PairTableSlot *oldSlots = table->slots;

	table->count = 0;
	table->capac = newCapac;
	table->slots = calloc(table->capac, sizeof *table->slots);

	for (size_t i = 0; i < oldCapac; i++) {
		if (oldSlots[i].occupied) {
			tin_insert_pair(table, oldSlots[i].elemLow, oldSlots[i].elemHigh, oldSlots[i].payload);
		}
	}

	free(oldSlots);
}

void
tin_order_pair(size_t *elemLow, size_t *elemHigh)
{
	if (*elemLow > *elemHigh) {
		size_t elemTmp = *elemLow;
		*elemLow = *elemHigh;
		*elemHigh = elemTmp;
	}
}

uint32_t
tin_hash_pair(size_t elemLow, size_t elemHigh)
{
	// TODO better hash function
	return ((elemLow ^ (elemHigh << 16)) * 33) >> 4;
}

size_t
tin_fold_hash(size_t capac, uint32_t hash)
{
	// TODO xor folding or similar
	uint32_t mask = capac - 1;
	return hash & mask;
}

size_t
tin_pairtable_index(Tin_PairTable *table, size_t elemLow, size_t elemHigh)
{
	uint32_t hash = tin_hash_pair(elemLow, elemHigh);
	size_t index = tin_fold_hash(table->capac, hash);
	for (;;) {
		if (!table->slots[index].occupied) {
			return index;
		}
		if (table->slots[index].hash == hash
		 && table->slots[index].elemLow == elemLow
		 && table->slots[index].elemHigh == elemHigh) {
			return index;
		}
		index++;
		if (index == table->capac) {
			index = 0;
		}
	}
}

bool
tin_find_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void **payloadOut)
{
	tin_order_pair(&elemA, &elemB);
	size_t index = tin_pairtable_index(table, elemA, elemB);
	if (!table->slots[index].occupied) {
		return false;
	}
	*payloadOut = table->slots[index].payload;
	return true;
}

void
tin_insert_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void *payload)
{
	if (3 * table->count >= 2 * table->capac) {
		tin_resize_pairtable(table, 2 * table->capac);
	}
	tin_order_pair(&elemA, &elemB);
	size_t index = tin_pairtable_index(table, elemA, elemB);
	if (!table->slots[index].occupied) {
		table->count++;
	}
	table->slots[index] = (Tin_PairTableSlot) {
		elemA,
		elemB,
		payload,
		tin_hash_pair(elemA, elemB),
		true,
	};
}

void
tin_delete_pair(Tin_PairTable *table, size_t elemA, size_t elemB)
{
	tin_order_pair(&elemA, &elemB);
	size_t index = tin_pairtable_index(table, elemA, elemB);
	if (!table->slots[index].occupied) {
		return;
	}
	table->count--;

	for (;;) {
		size_t next = index + 1;
		if (next == table->capac) {
			next = 0;
		}
		if (!table->slots[next].occupied) {
			break;
		}
		size_t nextDesired = tin_fold_hash(table->capac, table->slots[next].hash);
		if (next == nextDesired) {
			break;
		}
		table->slots[index] = table->slots[next];
		index = next;
	}

	table->slots[index].occupied = false;
}

/* === Island Detection & Freezing === :island: */

Tin_Body *
tin_island_find(Tin_Body *body)
{
	Tin_Body *island = body;
	while (island->island != island) {
		island = island->island;
	}
	if (body->island != island) {
		body->island = island;
	}
	return island;
}

void
tin_island_union(Tin_Body *bodyA, Tin_Body *bodyB)
{
	bodyA = tin_island_find(bodyA);
	bodyB = tin_island_find(bodyB);
	bodyB->island = bodyA;
}

void
tin_build_islands(Tin_Scene *scene, const Tin_Collision *collisions, size_t numCollisions)
{
	{
		TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
			body->island = body;
			body->islandStable = true;
		}
	}

	for (size_t c = 0; c < numCollisions; c++) {
		if (collisions[c].bodyA->invMass != 0.0 && collisions[c].bodyB->invMass != 0.0) {
			tin_island_union(collisions[c].bodyA, collisions[c].bodyB);
		}
	}

	{
		TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
			if (body->restCounter >= 5 || body->invMass == 0.0) continue;
			Tin_Body *island = tin_island_find(body);
			island->islandStable = false;
		}
	}
}

/* === Scenes / Worlds === :scene: */

Tin_Arbiter *
tin_find_arbiter(Tin_Scene *scene, Tin_Body *body1, Tin_Body *body2)
{
	void *payload;
	if (tin_find_pair(&scene->arbiters, (uintptr_t)body1, (uintptr_t)body2, &payload)) {
		return payload;
	}
	Tin_Arbiter *arbiter = tin_add_arbiter(scene);
	arbiter->body1 = body1;
	arbiter->body2 = body2;
	tin_insert_pair(&scene->arbiters, (uintptr_t)body1, (uintptr_t)body2, arbiter);
	return arbiter;
}

void
tin_scene_update(Tin_Scene *scene)
{
	TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
		Tin_Scalar *rotation = body->transform.rotation;

		/* Compute product = rotation * inertia_local^-1 */
		Tin_Scalar factor = body->invMass / (body->transform.scale * body->transform.scale);
		Tin_Scalar product[3*3];
		for (int j = 0; j < 3; j++) {
			Tin_Scalar diagonalEntry = factor * body->shape->invInertia.c[j];
			product[3*j+0] = rotation[3*j+0] * diagonalEntry;
			product[3*j+1] = rotation[3*j+1] * diagonalEntry;
			product[3*j+2] = rotation[3*j+2] * diagonalEntry;
		}

		/* Compute inertia_global^-1 = rotation * inertia_local^-1 * rotation^T */
		tin_m3_times_m3_transposed(body->invInertia, product, rotation);
	}

	for (size_t s = 0; s < scene->arbiters.capac; s++) {
		if (!scene->arbiters.slots[s].occupied) continue;
		tin_arbiter_update(scene->arbiters.slots[s].payload);
	}
}

void
tin_check_collision(Tin_Scene *scene, Tin_Body *body1, Tin_Body *body2)
{
	Tin_Arbiter *arbiter = tin_find_arbiter(scene, body1, body2);
	if (!arbiter) return;

	int numContacts = tin_polytope_collide(
			&body1->shape->polytope, &body1->transform,
			&body2->shape->polytope, &body2->transform, arbiter->contacts);
	arbiter->numContacts = numContacts;

	if (numContacts) {
		printf("%d contacts: ", numContacts);
		for (int c = 0; c < numContacts; c++) {
			const Tin_Contact *contact = &arbiter->contacts[c];
			printf(" (%f, %f, %f)", contact->position.c[0], contact->position.c[1], contact->position.c[2]);
		}
		printf("\n");
		/*printf("D %f\n", tin_dot_v3(contact.normal,
					tin_sub_v3(body2->transform.translation, body1->transform.translation)));*/
		//tin_arbiter_add_contact(arbiter, contact);
	}
}

void
tin_scene_prestep(Tin_Scene *scene, Tin_Collision *collisions, size_t numCollisions, Tin_Scalar invDt)
{
	for (size_t c = 0; c < numCollisions; c++) {
		Tin_Arbiter *arbiter = tin_find_arbiter(scene, collisions[c].bodyA, collisions[c].bodyB);
		tin_arbiter_prestep(arbiter, invDt);
		tin_arbiter_warm_start(arbiter);
	}
}

void
tin_scene_step(Tin_Scene *scene, Tin_Collision *collisions, size_t numCollisions, Tin_Scalar invDt)
{
	for (size_t c = 0; c < numCollisions; c++) {
		Tin_Arbiter *arbiter = tin_find_arbiter(scene, collisions[c].bodyA, collisions[c].bodyB);
		tin_arbiter_apply_separation(arbiter, invDt);
	}
	for (size_t c = 0; c < numCollisions; c++) {
		Tin_Arbiter *arbiter = tin_find_arbiter(scene, collisions[c].bodyA, collisions[c].bodyB);
		tin_arbiter_apply_friction(arbiter, invDt);
	}
	TIN_FOR_EACH(joint, scene->joints, Tin_Joint, node) {
		tin_joint_apply_impulse(joint, invDt);
	}
}

void
tin_integrate(Tin_Scene *scene, Tin_Scalar dt)
{
	TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
		bool stable = true;

		if (tin_dot_v3(body->velocity, body->velocity) > 10000.0f * TIN_EPSILON) {
			body->transform.translation = tin_saxpy_v3(dt, body->velocity, body->transform.translation);
			stable = false;
		}
		
		Tin_Scalar angle = sqrt(tin_dot_v3(body->angularVelocity, body->angularVelocity));
		
		if (fabs(angle) > 0.01) {
			Tin_Scalar change[3*3];
			tin_axis_angle_to_matrix(tin_normalize_v3(body->angularVelocity), angle * dt, change);

			Tin_Scalar copyOfRotation[3*3];
			memcpy(copyOfRotation, body->transform.rotation, sizeof copyOfRotation);
			
			tin_m3_times_m3(body->transform.rotation, change, copyOfRotation);

			stable = false;
		}

		if (!tin_island_find(body)->islandStable) {
			body->velocity = tin_saxpy_v3(dt, (Tin_Vec3) {{ 0.0f, -2.0f, 0.0f }}, body->velocity);
		}

		if (stable) {
			if (body->restCounter < 5) {
				body->restCounter++;
			}
		} else {
			body->restCounter = 0;
		}
	}
}

Tin_Collision *
tin_broadphase(Tin_Scene *scene, size_t *count_out)
{
	size_t capac = 0, count = 0;
	Tin_Collision *collisions = NULL;

	TIN_FOR_EACH(body1, scene->bodies, Tin_Body, node) {
		TIN_FOR_RANGE(body2, body1->node, scene->bodies, Tin_Body, node) {
			Tin_Vec3 diff = tin_sub_v3(body2->transform.translation, body1->transform.translation);
			Tin_Scalar radii =
				body2->shape->radius * body2->transform.scale +
				body1->shape->radius * body1->transform.scale;
			if (tin_dot_v3(diff, diff) <= radii * radii) {
				if (count == capac) {
					capac = capac ? 2*capac : 16;
					collisions = realloc(collisions, capac * sizeof *collisions);
				}
				collisions[count++] = (Tin_Collision){ body1, body2 };
			}
		}
	}

	*count_out = count;
	return collisions;
}

void
tin_narrowphase(Tin_Scene *scene, Tin_Collision *collisions, size_t num_collisions)
{
	for (size_t i = 0; i < num_collisions; i++) {
		tin_check_collision(scene, collisions[i].bodyA, collisions[i].bodyB);
	}
}

void
tin_simulate(Tin_Scene *scene, Tin_Scalar dt, double (*gettime)(), double timings[6])
{
	double startTime, stopTime;
	if (timings) {
		for (int t = 0; t < 6; t++) {
			timings[t] = 0.0;
		}
	}

	Tin_Scalar invDt = 1.0f / dt;

	startTime = gettime ? gettime() : 0.0;
	tin_scene_update(scene);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[0] += stopTime - startTime;
	
	startTime = gettime ? gettime() : 0.0;
	size_t num_collisions;
	Tin_Collision *collisions = tin_broadphase(scene, &num_collisions);
	tin_build_islands(scene, collisions, num_collisions);
	size_t old_num_collisions = num_collisions;
	num_collisions = 0;
	for (size_t c = 0; c < old_num_collisions; c++) {
		if (!tin_island_find(collisions[c].bodyA)->islandStable || !tin_island_find(collisions[c].bodyB)->islandStable) {
			collisions[num_collisions++] = collisions[c];
		}
	}
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[1] += stopTime - startTime;

	startTime = gettime ? gettime() : 0.0;
	tin_narrowphase(scene, collisions, num_collisions);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[2] += stopTime - startTime;
	
	startTime = gettime ? gettime() : 0.0;
	tin_scene_prestep(scene, collisions, num_collisions, invDt);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[3] += stopTime - startTime;

	startTime = gettime ? gettime() : 0.0;
	for (int iter = 0; iter < 16; iter++) {
		tin_scene_step(scene, collisions, num_collisions, invDt);
	}
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[4] += stopTime - startTime;
	
	startTime = gettime ? gettime() : 0.0;
	free(collisions);
	tin_integrate(scene, dt);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[5] += stopTime - startTime;
}

Tin_Body *
tin_add_body(Tin_Scene *scene, const Tin_Shape *shape, Tin_Scalar invMass)
{
	Tin_Body *body = scene->bodyAllocator.alloc(scene->bodyAllocator.userPointer);
	memset(body, 0, sizeof *body);
	body->transform = (Tin_Transform) {
		{ 1, 0, 0, 0, 1, 0, 0, 0, 1 },
		TIN_VEC3(0, 0, 0),
		1.0,
	};
	body->shape = shape;
	body->invMass = invMass;
	TIN_LIST_LINK(*scene->bodies.prev, body->node);
	TIN_LIST_LINK(body->node, scene->bodies);
	return body;
}

Tin_Arbiter *
tin_add_arbiter(Tin_Scene *scene)
{
	Tin_Arbiter *arbiter = scene->arbiterAllocator.alloc(scene->arbiterAllocator.userPointer);
	memset(arbiter, 0, sizeof *arbiter);
	return arbiter;
}

Tin_Joint *
tin_add_joint(Tin_Scene *scene)
{
	Tin_Joint *joint = scene->jointAllocator.alloc(scene->jointAllocator.userPointer);
	memset(joint, 0, sizeof *joint);
	TIN_LIST_LINK(*scene->joints.prev, joint->node);
	TIN_LIST_LINK(joint->node, scene->joints);
	return joint;
}

/* === Broadphase === :broad: */

#if 0
struct Endpoint {
	size_t bodyIdx;
	Tin_Scalar value;
	bool isMin;
};

struct Broadphase {
	size_t numBodies;
	struct Endpoint *axes[3];
};

void
tin_create_broadphase(struct Broadphase *broad, size_t numBodies)
{
	broad->numBodies = numBodies;
	broad->touching = calloc(numBodies, sizeof *broad->touching);
	broad->axes[0] = calloc(2 * numBodies, sizeof *broad->axes[0]);
	broad->axes[1] = calloc(2 * numBodies, sizeof *broad->axes[0]);
	broad->axes[2] = calloc(2 * numBodies, sizeof *broad->axes[0]);
}

void
tin_destroy_broadphase(struct Broadphase *broad)
{
	for (size_t b = 0; b < broad->numBodies; b++) {
		struct BodyList *list = broad->touching[b];
		while (list) {
			struct BodyList *next = list->next;
			free(list);
			list = next;
		}
	}
	free(broad->touching);
	free(broad->axes[0]);
	free(broad->axes[1]);
	free(broad->axes[2]);
}

void
tin_update_broadphase(struct Broadphase *broad)
{
	(void) broad;
	for () {
		Tin_Body *body = tin_body_pointer(scene, index);
	}
}
#endif

#if 0
void
tin_sort_axis(struct Broadphase *broad, int axisNum)
{
	struct Endpoint *axis = broad->axes[axisNum];

	for (size_t j, i = 1; i < 2 * broad->numBodies; i++) {
		struct Endpoint point = axis[i];
		for (j = i; j > 0 && axis[j-1].value > point.value; j--) {
			axis[j] = axis[j-1];
		}
		axis[j] = point;
	}
}

void
tin_scan_axis(struct Broadphase *broad, int axisNum)
{
	struct Endpoint *axis = broad->axes[axisNum];
	size_t numEndpoints = 2 * broad->numBodies;

	size_t *active = calloc(numEndpoints, sizeof *active);
	size_t numActive = 0;

	for (size_t i = 0; i < numEndpoints; i++) {
		if (axis[i].isMin) {
			for (size_t a = 0; a < numActive; a++) {
				size_t bodyA = active[a];
				size_t bodyB = axis[i].bodyIdx;

				struct BodyList *listNodeA = calloc(1, sizeof *listNodeA);
				listNodeA->bodyIdx = bodyB;
				listNodeA->next = broad->touching[bodyA];
				broad->touching[bodyA] = listNodeA;

				struct BodyList *listNodeB = calloc(1, sizeof *listNodeB);
				listNodeB->bodyIdx = bodyA;
				listNodeB->next = broad->touching[bodyB];
				broad->touching[bodyB] = listNodeB;
			}
			active[numActive++] = axis[i].bodyIdx;
		} else {
			for (size_t a = 0; a < numActive; a++) {
				if (active[a] == axis[i].bodyIdx) {
					active[a] = active[--numActive];
					break;
				}
			}
		}
	}

	free(active);
}
#endif

