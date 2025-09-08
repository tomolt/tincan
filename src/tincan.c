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

#ifdef _MSC_VER // compiling with MSVC
#	if defined(_M_AMD64) || defined(_M_X64) || _M_IX86_FP == 2
#		define TIN_HAS_SSE2 1
#		include <intrin.h>
#	endif
#else // any other compiler; Assume compatibility with GCC
#	if defined(__SSE2__)
#		define TIN_HAS_SSE2 1
#		include <emmintrin.h>
#	endif
#endif

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* === Vector Math === :vec: */

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
	// Clip planes point outwards. Polytope faces need to be wound counter-clockwise.
	int newCount = countB;
	memcpy(newPoints, pointsB, newCount * sizeof *newPoints);
	int j = countA - 1;
	for (int i = 0; i < countA; i++) {
		if (newCount == 0) break;
		Tin_Vec3 edge = tin_sub_v3(pointsA[i], pointsA[j]);
		Tin_Vec3 normal = tin_cross_v3(edge, perpendicular);
		normal = tin_normalize_v3(normal);
		Tin_Scalar base = tin_dot_v3(normal, pointsA[i]);
		Tin_Vec3 buffer[32];
		newCount = tin_clip_manifold_against_plane(newPoints, newCount, normal, base, buffer);
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
	Tin_Arbiter *arbiter)
{
	Tin_Polysum ps = { pa, ta, pb, tb };

	Tin_Ray r;
	r.origin = tin_sub_v3(ta->translation, tb->translation);
	r.dir    = tin_neg_v3(r.origin);
	Tin_Scalar norm = sqrt(tin_dot_v3(r.dir, r.dir));
	if (norm == 0.0) {
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
	if (tin_dot_v3(p.normal, p.a) <= 0.0) {
		return 0;
	}

	for (int it = 0; it < 4; it++) {
		Tin_Ray nr;
		nr.dir    = p.normal;
		nr.origin = TIN_VEC3(0.0, 0.0, 0.0);
		Tin_Portal np;
		if (!tin_construct_portal(&ps, tin_polysum_support, &nr, &np)) {
			break;
		}
		tin_refine_portal(&ps, tin_polysum_support, &nr, &np);
		np.normal = tin_normalize_v3(np.normal);

		Tin_Scalar proj = tin_dot_v3(p.normal, np.normal);
		p = np; /* lol */
		r = nr;
		if (proj >= 0.99) break;
	}

	int faceA = tin_incident_face(pa, tin_bwtrf_dir(ta, p.normal));
	Tin_Vec3 pointsA[32];
	int countA = 0;
	for (int i = pa->faceOffsets[faceA]; i < pa->faceOffsets[faceA+1]; i++) {
		int idx = pa->faceIndices[i];
		pointsA[countA++] = tin_fwtrf_point(ta, pa->vertices[idx]);
	}

	Tin_Vec3 refNormal = tin_normalize_v3(tin_fwtrf_dir(ta, pa->faceNormals[faceA]));
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

	arbiter->normal = refNormal;
	arbiter->face1 = faceA;
	arbiter->face2 = faceB;
	memset(arbiter->contacts, 0, count * sizeof *arbiter->contacts);
	for (int i = 0; i < count; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];
		Tin_Scalar projDist = tin_dot_v3(refNormal, manifold[i]) - refBase;
		Tin_Vec3 projPoint = tin_saxpy_v3(-projDist, refNormal, manifold[i]);
		contact->posFrom1 = tin_sub_v3(projPoint, ta->translation);
		contact->posFrom2 = tin_sub_v3(manifold[i], tb->translation);
		Tin_Vec3 p1 = projPoint;
		Tin_Vec3 p2 = manifold[i];
		contact->separation = tin_dot_v3(refNormal, tin_sub_v3(p2, p1));
	}

	for (int n = count; n > 1; n--) {
		for (int i = 0; i < n - 1; i++) {
			if (arbiter->contacts[i].separation >= 0.0 && !(arbiter->contacts[i+1].separation >= 0.0)) {
				Tin_Contact temp = arbiter->contacts[i];
				arbiter->contacts[i] = arbiter->contacts[i+1];
				arbiter->contacts[i+1] = temp;
			}
		}
	}

	int numPenetrating = 0;
	for (int i = 0; i < count; i++) {
		if (!(arbiter->contacts[i].separation >= 0.0)) {
			numPenetrating++;
		}
	}
	arbiter->numPenetrating = numPenetrating;

	return count;
}

Tin_Scalar
tin_dot_v12(Tin_Scalar a[12], Tin_Scalar b[12])
{
#if TIN_HAS_SSE2
	float f[4];
	__m128 x0, x1, x2;
	x0 = _mm_mul_ps(_mm_loadu_ps(&a[0]), _mm_loadu_ps(&b[0]));
	x1 = _mm_mul_ps(_mm_loadu_ps(&a[4]), _mm_loadu_ps(&b[4]));
	x2 = _mm_mul_ps(_mm_loadu_ps(&a[8]), _mm_loadu_ps(&b[8]));
	x0 = _mm_add_ps(x0, x1);
	x0 = _mm_add_ps(x0, x2);
	_mm_store_ps(f, x0);
	return (f[0] + f[1]) + (f[2] + f[3]);
#else
	Tin_Scalar sum = 0.0;
	for (int i = 0; i < 12; i++) {
		sum += a[i] * b[i];
	}
	return sum;
#endif
}

Tin_Scalar
tin_dot_array(Tin_Scalar a[], Tin_Scalar b[], int length)
{
	Tin_Scalar sum = 0.0;
	for (int i = 0; i < length; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

Tin_Scalar
tin_effective_mass(const Tin_Body *body1, const Tin_Body *body2, Tin_Scalar jacobian[12])
{
	Tin_Scalar invMassJacobian[12];
	Tin_Vec3 rotMass1 = tin_m3_times_v3(body1->invInertia, TIN_VEC3(jacobian[ 2], jacobian[ 6], jacobian[10]));
	Tin_Vec3 rotMass2 = tin_m3_times_v3(body2->invInertia, TIN_VEC3(jacobian[ 3], jacobian[ 7], jacobian[11]));
	invMassJacobian[ 0] = body1->invMass * jacobian[0];
	invMassJacobian[ 1] = body2->invMass * jacobian[1];
	invMassJacobian[ 2] = rotMass1.c[0];
	invMassJacobian[ 3] = rotMass2.c[0];
	invMassJacobian[ 4] = body1->invMass * jacobian[4];
	invMassJacobian[ 5] = body2->invMass * jacobian[5];
	invMassJacobian[ 6] = rotMass1.c[1];
	invMassJacobian[ 7] = rotMass2.c[1];
	invMassJacobian[ 8] = body1->invMass * jacobian[8];
	invMassJacobian[ 9] = body2->invMass * jacobian[9];
	invMassJacobian[10] = rotMass1.c[2];
	invMassJacobian[11] = rotMass2.c[2];

	Tin_Scalar effectiveMass = 1.0 / tin_dot_v12(jacobian, invMassJacobian);
	return effectiveMass;
}

// OUTDATED
void
tin_apply_impulse(Tin_Body *body1, Tin_Body *body2, Tin_Scalar jacobian[12], Tin_Scalar magnitude)
{
	abort();
	body1->velocity = tin_saxpy_v3(magnitude * body1->invMass, TIN_VEC3(jacobian[0], jacobian[1], jacobian[2]), body1->velocity);
	body1->angularVelocity = tin_saxpy_v3(magnitude, tin_m3_times_v3(body1->invInertia, TIN_VEC3(jacobian[3], jacobian[4], jacobian[5])), body1->angularVelocity);
	body2->velocity = tin_saxpy_v3(magnitude * body2->invMass, TIN_VEC3(jacobian[6], jacobian[7], jacobian[8]), body2->velocity);
	body2->angularVelocity = tin_saxpy_v3(magnitude, tin_m3_times_v3(body2->invInertia, TIN_VEC3(jacobian[9], jacobian[10], jacobian[11])), body2->angularVelocity);
}

Tin_Scalar
tin_enforce_jacobian_fast(Tin_Scalar (*velocities)[6], int body1Idx, int body2Idx, Tin_Scalar body1InvMass, Tin_Scalar body2InvMass, Tin_Scalar jacobian[12],
	Tin_Vec3 angularImpulse1, Tin_Vec3 angularImpulse2,
	Tin_Scalar effectiveMass, Tin_Scalar bias,
	Tin_Scalar magnitudeAccum, Tin_Scalar minMagnitude, Tin_Scalar maxMagnitude)
{
	/* Solve for magnitude */
	Tin_Scalar velocity[12];
	velocity[ 0] = velocities[body1Idx][0];
	velocity[ 1] = velocities[body2Idx][0];
	velocity[ 2] = velocities[body1Idx][1];
	velocity[ 3] = velocities[body2Idx][1];
	velocity[ 4] = velocities[body1Idx][2];
	velocity[ 5] = velocities[body2Idx][2];
	velocity[ 6] = velocities[body1Idx][3];
	velocity[ 7] = velocities[body2Idx][3];
	velocity[ 8] = velocities[body1Idx][4];
	velocity[ 9] = velocities[body2Idx][4];
	velocity[10] = velocities[body1Idx][5];
	velocity[11] = velocities[body2Idx][5];
	Tin_Scalar magnitude = bias;
	magnitude -= tin_dot_v12(jacobian, velocity);
	magnitude *= effectiveMass;

	/* Clamp magnitude to fulfill inequality */
	Tin_Scalar prevAccum = magnitudeAccum;
	magnitudeAccum += magnitude;
	magnitudeAccum = MAX(magnitudeAccum, minMagnitude);
	magnitudeAccum = MIN(magnitudeAccum, maxMagnitude);
	magnitude = magnitudeAccum - prevAccum;

#if TIN_HAS_SSE2
	__m128 vx = _mm_set_ps(velocities[body2Idx][1], velocities[body2Idx][0], velocities[body1Idx][1], velocities[body1Idx][0]);
	__m128 vy = _mm_set_ps(velocities[body2Idx][3], velocities[body2Idx][2], velocities[body1Idx][3], velocities[body1Idx][2]);
	__m128 vz = _mm_set_ps(velocities[body2Idx][5], velocities[body2Idx][4], velocities[body1Idx][5], velocities[body1Idx][4]);
	__m128 factor = _mm_set1_ps(magnitude) * _mm_set_ps(1.0, body2InvMass, 1.0, body1InvMass);
	vx = _mm_add_ps(vx, _mm_mul_ps(factor, _mm_set_ps(angularImpulse2.c[0], jacobian[1], angularImpulse1.c[0], jacobian[0])));
	vy = _mm_add_ps(vy, _mm_mul_ps(factor, _mm_set_ps(angularImpulse2.c[1], jacobian[5], angularImpulse1.c[1], jacobian[4])));
	vz = _mm_add_ps(vz, _mm_mul_ps(factor, _mm_set_ps(angularImpulse2.c[2], jacobian[9], angularImpulse1.c[2], jacobian[8])));
	float f[12];
	_mm_store_ps(&f[0], vx);
	_mm_store_ps(&f[4], vy);
	_mm_store_ps(&f[8], vz);
	velocities[body1Idx][0] = f[0];
	velocities[body1Idx][1] = f[1];
	velocities[body2Idx][0] = f[2];
	velocities[body2Idx][1] = f[3];
	velocities[body1Idx][2] = f[4];
	velocities[body1Idx][3] = f[5];
	velocities[body2Idx][2] = f[6];
	velocities[body2Idx][3] = f[7];
	velocities[body1Idx][4] = f[8];
	velocities[body1Idx][5] = f[9];
	velocities[body2Idx][4] = f[10];
	velocities[body2Idx][5] = f[11];
#else
	// TODO Non-SIMD tin_enforce_jacobian_fast()
#error "Non-SIMD tin_enforce_jacobian_fast() is not implemented yet"
#endif

	return magnitudeAccum;
}

void
tin_enforce_jacobian(Tin_Body *body1, Tin_Body *body2, Tin_Scalar jacobian[12],
	Tin_Scalar effectiveMass, Tin_Scalar bias,
	Tin_Scalar *magnitudeAccum, Tin_Scalar minMagnitude, Tin_Scalar maxMagnitude)
{
	abort();

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
tin_jacobian_along_axis(Tin_Scalar jacobian[12], Tin_Vec3 axis, Tin_Vec3 r1, Tin_Vec3 r2)
{
	Tin_Vec3 r1Xaxis = tin_cross_v3(r1, axis);
	Tin_Vec3 r2Xaxis = tin_cross_v3(r2, axis);
	// SIMD friendly layout: V1X V2X W1X W2X V1Y V2Y W1Y W2Y V1Z V2Z W1Z W2Z
	jacobian[ 0] = -axis.c[0];
	jacobian[ 1] =  axis.c[0];
	jacobian[ 2] = -r1Xaxis.c[0];
	jacobian[ 3] =  r2Xaxis.c[0];
	jacobian[ 4] = -axis.c[1];
	jacobian[ 5] =  axis.c[1];
	jacobian[ 6] = -r1Xaxis.c[1];
	jacobian[ 7] =  r2Xaxis.c[1];
	jacobian[ 8] = -axis.c[2];
	jacobian[ 9] =  axis.c[2];
	jacobian[10] = -r1Xaxis.c[2];
	jacobian[11] =  r2Xaxis.c[2];
}

void
tin_arbiter_prestep(Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6], Tin_Scalar invDt)
{
	const Tin_Scalar allowedPenetration = 0.01;
	const Tin_Scalar biasFactor = 0.1;

	Tin_Scalar *vel1 = velocities[arbiter->body1Idx];
	Tin_Scalar *vel2 = velocities[arbiter->body2Idx];
	Tin_Vec3 relVel = tin_sub_v3(TIN_VEC3(vel1[0], vel1[2], vel1[4]), TIN_VEC3(vel2[0], vel2[2], vel2[4]));
	Tin_Vec3 frictionDir = tin_gram_schmidt(arbiter->normal, relVel);
	if (tin_dot_v3(frictionDir, frictionDir) < 0.001) {
		frictionDir = tin_gram_schmidt(arbiter->normal, TIN_VEC3(1.0, 0.0, 0.0));
		if (tin_dot_v3(frictionDir, frictionDir) == 0.0) {
			frictionDir = tin_gram_schmidt(arbiter->normal, TIN_VEC3(0.0, 1.0, 0.0));
		}
	}
	frictionDir = tin_normalize_v3(frictionDir);
	Tin_Vec3 orthoDir = tin_normalize_v3(tin_cross_v3(arbiter->normal, frictionDir));

	arbiter->frictionDir = frictionDir;
	arbiter->orthoDir = orthoDir;

	arbiter->body1InvMass = arbiter->body1->invMass;
	arbiter->body2InvMass = arbiter->body2->invMass;

	for (int i = 0; i < arbiter->numContacts; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		/* Precompute jacobian, effectiveMass, and bias */
		tin_jacobian_along_axis(contact->jacobian, arbiter->normal, contact->posFrom1, contact->posFrom2);
		contact->effectiveMass[0] = tin_effective_mass(arbiter->body1, arbiter->body2, contact->jacobian);

		contact->bias = -biasFactor * invDt * MIN(0.0, contact->separation + allowedPenetration);

		tin_jacobian_along_axis(contact->tangentJacobian, frictionDir, contact->posFrom1, contact->posFrom2);
		contact->effectiveMass[1] = tin_effective_mass(arbiter->body1, arbiter->body2, contact->tangentJacobian);
		contact->tangentAngularImpulse1 = tin_m3_times_v3(arbiter->body1->invInertia,
			TIN_VEC3(contact->tangentJacobian[2], contact->tangentJacobian[6], contact->tangentJacobian[10]));
		contact->tangentAngularImpulse2 = tin_m3_times_v3(arbiter->body2->invInertia,
			TIN_VEC3(contact->tangentJacobian[3], contact->tangentJacobian[7], contact->tangentJacobian[11]));


		tin_jacobian_along_axis(contact->bitangentJacobian, orthoDir, contact->posFrom1, contact->posFrom2);
		contact->effectiveMass[2] = tin_effective_mass(arbiter->body1, arbiter->body2, contact->bitangentJacobian);
		contact->bitangentAngularImpulse1 = tin_m3_times_v3(arbiter->body1->invInertia,
			TIN_VEC3(contact->bitangentJacobian[2], contact->bitangentJacobian[6], contact->bitangentJacobian[10]));
		contact->bitangentAngularImpulse2 = tin_m3_times_v3(arbiter->body2->invInertia,
			TIN_VEC3(contact->bitangentJacobian[3], contact->bitangentJacobian[7], contact->bitangentJacobian[11]));


		contact->normalAngularImpulse1 = tin_m3_times_v3(arbiter->body1->invInertia,
			TIN_VEC3(contact->jacobian[2], contact->jacobian[6], contact->jacobian[10]));
		contact->normalAngularImpulse2 = tin_m3_times_v3(arbiter->body2->invInertia,
			TIN_VEC3(contact->jacobian[3], contact->jacobian[7], contact->jacobian[11]));

		contact->ineqAccum = 0.0;
		contact->tangentAccum = 0.0;
		contact->bitangentAccum = 0.0;
	}
}

void
tin_arbiter_warm_start(Tin_Arbiter *arbiter, const Tin_Arbiter *oldArbiter)
{
	const Tin_Scalar maxDistanceSq = 0.01;
	if (!(arbiter->face1 == oldArbiter->face1 && arbiter->face2 == oldArbiter->face2)) {
		return;
	}
	for (int i = 0; i < arbiter->numPenetrating; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		for (int j = 0; j < oldArbiter->numContacts; j++) {
			const Tin_Contact *oldContact = &oldArbiter->contacts[j];
			Tin_Vec3 diff1 = tin_sub_v3(contact->posFrom1, oldContact->posFrom1);
			Tin_Vec3 diff2 = tin_sub_v3(contact->posFrom2, oldContact->posFrom2);
			Tin_Scalar dist1Sq = tin_dot_v3(diff1, diff1);
			Tin_Scalar dist2Sq = tin_dot_v3(diff2, diff2);
			if (dist1Sq <= maxDistanceSq && dist2Sq <= maxDistanceSq) {
				contact->ineqAccum = oldContact->ineqAccum;
				contact->tangentAccum = oldContact->tangentAccum;
				contact->bitangentAccum = oldContact->bitangentAccum;
				tin_apply_impulse(arbiter->body1, arbiter->body2, contact->jacobian, contact->ineqAccum);
				tin_apply_impulse(arbiter->body1, arbiter->body2, contact->tangentJacobian, contact->tangentAccum);
				tin_apply_impulse(arbiter->body1, arbiter->body2, contact->bitangentJacobian, contact->bitangentAccum);
				break;
			}
		}
	}
}

void
tin_arbiter_apply_separation(Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6], Tin_Scalar invDt)
{
	// TODO
	(void)invDt;

	for (int idx = 0; idx < arbiter->numPenetrating; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];

		contact->ineqAccum = tin_enforce_jacobian_fast(velocities, arbiter->body1Idx, arbiter->body2Idx, arbiter->body1InvMass, arbiter->body2InvMass, contact->jacobian, contact->normalAngularImpulse1, contact->normalAngularImpulse2, contact->effectiveMass[0], contact->bias, contact->ineqAccum, 0.0, INFINITY);
	}
}

void
tin_arbiter_apply_friction(Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6], Tin_Scalar invDt)
{
	// TODO
	(void)invDt;

	const Tin_Scalar friction = 0.4;

	for (int idx = 0; idx < arbiter->numPenetrating; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];

		contact->tangentAccum = tin_enforce_jacobian_fast(velocities, arbiter->body1Idx, arbiter->body2Idx, arbiter->body1InvMass, arbiter->body2InvMass, contact->tangentJacobian, contact->tangentAngularImpulse1, contact->tangentAngularImpulse2, contact->effectiveMass[1],
			0.0, contact->tangentAccum, -contact->ineqAccum * friction, contact->ineqAccum * friction);

		contact->bitangentAccum = tin_enforce_jacobian_fast(velocities, arbiter->body1Idx, arbiter->body2Idx, arbiter->body1InvMass, arbiter->body2InvMass, contact->bitangentJacobian, contact->bitangentAngularImpulse1, contact->bitangentAngularImpulse2, contact->effectiveMass[2],
			0.0, contact->bitangentAccum, -contact->ineqAccum * friction, contact->ineqAccum * friction);
	}
}

void
tin_joint_apply_impulse(Tin_Joint *joint, Tin_Scalar invDt)
{
	const Tin_Scalar biasFactor = 0.1;

	Tin_Scalar bias;
	Tin_Scalar separation;
	Tin_Scalar jacobian[12];
	Tin_Scalar effectiveMass;

	abort();

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
	tin_enforce_jacobian(joint->body1, joint->body2, jacobian, effectiveMass, bias, NULL, 0.0, 0.0);

	separation = p2.c[1] - p1.c[1];
	bias = -biasFactor * invDt * separation;
	tin_jacobian_along_axis(jacobian, (Tin_Vec3){{0.0f, 1.0f, 0.0f}}, r1, r2);
	effectiveMass = tin_effective_mass(joint->body1, joint->body2, jacobian);
	tin_enforce_jacobian(joint->body1, joint->body2, jacobian, effectiveMass, bias, NULL, 0.0, 0.0);

	separation = p2.c[2] - p1.c[2];
	bias = -biasFactor * invDt * separation;
	tin_jacobian_along_axis(jacobian, (Tin_Vec3){{0.0f, 0.0f, 1.0f}}, r1, r2);
	effectiveMass = tin_effective_mass(joint->body1, joint->body2, jacobian);
	tin_enforce_jacobian(joint->body1, joint->body2, jacobian, effectiveMass, bias, NULL, 0.0, 0.0);
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
tin_island_union(Tin_Body *body1, Tin_Body *body2)
{
	body1 = tin_island_find(body1);
	body2 = tin_island_find(body2);
	body2->island = body1;
}

void
tin_build_islands(Tin_Scene *scene)
{
	{
		TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
			body->island = body;
			body->islandStable = true;
		}
	}

	for (size_t i = 0; i < scene->numArbiters; i++) {
		Tin_Arbiter *arbiter = &scene->arbiters[i];
		if (arbiter->body1->invMass != 0.0 && arbiter->body2->invMass != 0.0) {
			tin_island_union(arbiter->body1, arbiter->body2);
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
}

void
tin_scene_prestep(Tin_Scene *scene, Tin_Scalar (*velocities)[6], Tin_Scalar invDt)
{
	for (size_t i = 0; i < scene->numArbiters; i++) {
		tin_arbiter_prestep(&scene->arbiters[i], velocities, invDt);
#if 0
		void *payload;
		if (tin_find_pair(&scene->contactCache, (uintptr_t)scene->arbiters[i].body1, (uintptr_t)scene->arbiters[i].body2, &payload)) {
			tin_arbiter_warm_start(&scene->arbiters[i], payload);
		}
#endif
	}
}

void
tin_scene_step(Tin_Scene *scene, Tin_Scalar (*velocities)[6], Tin_Scalar invDt)
{
	for (size_t i = 0; i < scene->numArbiters; i++) {
		tin_arbiter_apply_separation(&scene->arbiters[i], velocities, invDt);
	}
	for (size_t i = 0; i < scene->numArbiters; i++) {
		tin_arbiter_apply_friction(&scene->arbiters[i], velocities, invDt);
	}
#if 0
	TIN_FOR_EACH(joint, scene->joints, Tin_Joint, node) {
		tin_joint_apply_impulse(joint, invDt);
	}
#endif
}

void
tin_integrate(Tin_Scene *scene, Tin_Scalar dt)
{
	TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
		bool stable = true;

		if (tin_dot_v3(body->velocity, body->velocity) > 1000.0f * TIN_EPSILON) {
			body->transform.translation = tin_saxpy_v3(dt, body->velocity, body->transform.translation);
			stable = false;
		} else {
			body->velocity = TIN_VEC3(0, 0, 0);
		}
		
		Tin_Scalar angle = sqrt(tin_dot_v3(body->angularVelocity, body->angularVelocity));
		
		if (fabs(angle) > 0.02) {
			Tin_Scalar change[3*3];
			tin_axis_angle_to_matrix(tin_normalize_v3(body->angularVelocity), angle * dt, change);

			Tin_Scalar copyOfRotation[3*3];
			memcpy(copyOfRotation, body->transform.rotation, sizeof copyOfRotation);
			
			tin_m3_times_m3(body->transform.rotation, change, copyOfRotation);

			stable = false;
		} else {
			body->angularVelocity = TIN_VEC3(0, 0, 0);
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

void
tin_broadphase(Tin_Scene *scene)
{
	size_t capac = scene->capArbiters, count = 0;
	Tin_Arbiter *arbiters = scene->arbiters;

	int body1Idx = 0;
	TIN_FOR_EACH(body1, scene->bodies, Tin_Body, node) {
		int body2Idx = body1Idx + 1;
		TIN_FOR_RANGE(body2, body1->node, scene->bodies, Tin_Body, node) {
			Tin_Vec3 diff = tin_sub_v3(body2->transform.translation, body1->transform.translation);
			Tin_Scalar radii =
				body2->shape->radius * body2->transform.scale +
				body1->shape->radius * body1->transform.scale;
			if (tin_dot_v3(diff, diff) <= radii * radii) {
				if (count == capac) {
					capac = capac ? 2*capac : 16;
					arbiters = realloc(arbiters, capac * sizeof *arbiters);
				}
				arbiters[count++] = (Tin_Arbiter){ .body1 = body1, .body2 = body2, .body1Idx = body1Idx, .body2Idx = body2Idx };
			}
			body2Idx++;
		}
		body1Idx++;
	}

	scene->arbiters = arbiters;
	scene->numArbiters = count;
	scene->capArbiters = capac;
}

void
tin_narrowphase(Tin_Scene *scene)
{
	for (size_t i = 0; i < scene->numArbiters; i++) {
		scene->arbiters[i].numContacts = tin_polytope_collide(
				&scene->arbiters[i].body1->shape->polytope, &scene->arbiters[i].body1->transform,
				&scene->arbiters[i].body2->shape->polytope, &scene->arbiters[i].body2->transform, &scene->arbiters[i]);
	}
	size_t numOut = 0;
	for (size_t i = 0; i < scene->numArbiters; i++) {
		if (scene->arbiters[i].numContacts > 0) {
			scene->arbiters[numOut++] = scene->arbiters[i];
		}
	}
	scene->numArbiters = numOut;
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
	tin_broadphase(scene);
	tin_build_islands(scene);
	/* Apply Gravity */
	size_t numBodies = 0;
	{
		TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
			if (!tin_island_find(body)->islandStable) {
				body->velocity = tin_saxpy_v3(dt, TIN_VEC3(0.0, -9.0, 0.0), body->velocity);
			}
			numBodies++;
		}
	}

	/* Filter out collisions within stable islands */
	size_t oldNumArbiters = scene->numArbiters;
	scene->numArbiters = 0;
	for (size_t i = 0; i < oldNumArbiters; i++) {
		if (tin_island_find(scene->arbiters[i].body1)->islandStable) {
			if (tin_island_find(scene->arbiters[i].body2)->islandStable) {
				continue;
			}
		}
		scene->arbiters[scene->numArbiters++] = scene->arbiters[i];
	}
	printf("#collisions = %zu\n", scene->numArbiters);

	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[1] += stopTime - startTime;

	startTime = gettime ? gettime() : 0.0;
	tin_narrowphase(scene);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[2] += stopTime - startTime;
	
	// We keep some padding at the end so we can more easily load values into SIMD registers.
	Tin_Scalar (*velocities)[6] = calloc((numBodies + 1), sizeof (Tin_Scalar[6]));
	{
		size_t b = 0;
		TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
			velocities[b][0] = body->velocity.c[0];
			velocities[b][1] = body->angularVelocity.c[0];
			velocities[b][2] = body->velocity.c[1];
			velocities[b][3] = body->angularVelocity.c[1];
			velocities[b][4] = body->velocity.c[2];
			velocities[b][5] = body->angularVelocity.c[2];
			b++;
		}
	}
	
	startTime = gettime ? gettime() : 0.0;
	tin_scene_prestep(scene, velocities, invDt);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[3] += stopTime - startTime;

	startTime = gettime ? gettime() : 0.0;
	for (int iter = 0; iter < 16; iter++) {
		tin_scene_step(scene, velocities, invDt);
	}
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[4] += stopTime - startTime;

	{
		size_t b = 0;
		TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
			body->velocity.c[0] = velocities[b][0];
			body->angularVelocity.c[0] = velocities[b][1];
			body->velocity.c[1] = velocities[b][2];
			body->angularVelocity.c[1] = velocities[b][3];
			body->velocity.c[2] = velocities[b][4];
			body->angularVelocity.c[2] = velocities[b][5];
			b++;
		}
	}
	free(velocities);

	startTime = gettime ? gettime() : 0.0;
	// Cache this frames contacts for the next frame
	/*
	for (size_t i = 0; i < numArbiters; i++) {
		Tin_CachedContact *cachedContact = calloc(1, sizeof *cachedContact);
		cachedContact->face1 = arbiters[i].face1;
		cachedContact->face2 = arbiters[i].face2;
		cachedContact->numPoints = arbiters[i].numContacts;
		for (int p = 0; p < arbiters[i].numContacts; p++) {
			Tin_CachedContactPoint *cachedPoint = &cachedContact->points[p];
			cachedPoint->posFrom1 = arbiters[i].contacts[p].posFrom1;
			cachedPoint->posFrom2 = arbiters[i].contacts[p].posFrom2;
			cachedPoint->magnitudeAccums[0] = arbiters[i].contacts[p].ineqAccum;
			cachedPoint->magnitudeAccums[1] = arbiters[i].contacts[p].tangentAccum;
			cachedPoint->magnitudeAccums[2] = arbiters[i].contacts[p].bitangentAccum;
		}
		tin_insert_pair(&scene->contactCache, (uintptr_t)arbiters[i].body1, (uintptr_t)arbiters[i].body2, cachedContact);
	}
	*/
	scene->numOldArbiters = scene->numArbiters;
	scene->capOldArbiters = scene->capArbiters;
	scene->oldArbiters = realloc(scene->oldArbiters, scene->capOldArbiters * sizeof *scene->oldArbiters);
	memcpy(scene->oldArbiters, scene->arbiters, scene->numOldArbiters * sizeof *scene->oldArbiters);
	for (size_t i = 0; i < scene->numOldArbiters; i++) {
		Tin_Arbiter *arbiter = &scene->oldArbiters[i];
		tin_insert_pair(&scene->contactCache,
			(uintptr_t)arbiter->body1, (uintptr_t)arbiter->body2,
			arbiter);
	}
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
				size_t body1 = active[a];
				size_t body2 = axis[i].bodyIdx;

				struct BodyList *listNodeA = calloc(1, sizeof *listNodeA);
				listNodeA->bodyIdx = body2;
				listNodeA->next = broad->touching[body1];
				broad->touching[body1] = listNodeA;

				struct BodyList *listNodeB = calloc(1, sizeof *listNodeB);
				listNodeB->bodyIdx = body1;
				listNodeB->next = broad->touching[body2];
				broad->touching[body2] = listNodeB;
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

