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
#		include <stdalign.h>
#		include <intrin.h>
#	endif
#else // any other compiler; Assume compatibility with GCC
#	if defined(__SSE2__)
#		define TIN_HAS_SSE2 1
#		include <stdalign.h>
#		include <emmintrin.h>
#	endif
#endif

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

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
		// A score of INFINITY is a distinct possibility (e1 = 0),
		// so we have to make sure that we still pick anybody in that
		// case.
		if (score <= bestScore) {
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
	Tin_Vec3 normal;
	if (!tin_intersect(pa, ta, pb, tb, &normal)) {
		return 0;
	}

	// Find a face on polytope A that aligns closely with the collision normal.
	int faceA = tin_incident_face(pa, tin_bwtrf_dir(ta, normal));
	Tin_Vec3 pointsA[32];
	int countA = 0;
	for (int i = pa->faceOffsets[faceA]; i < pa->faceOffsets[faceA+1]; i++) {
		int idx = pa->faceIndices[i];
		pointsA[countA++] = tin_fwtrf_point(ta, pa->vertices[idx]);
	}

	Tin_Vec3 refNormal = tin_normalize_v3(tin_fwtrf_dir(ta, pa->faceNormals[faceA]));
	Tin_Scalar refBase = tin_dot_v3(refNormal, pointsA[0]);

	// Find a face on polytope B that aligns closely with the (negated) collision normal.
	int faceB = tin_incident_face(pb, tin_bwtrf_dir(tb, tin_neg_v3(refNormal)));
	Tin_Vec3 pointsB[32];
	int countB = 0;
	for (int i = pb->faceOffsets[faceB]; i < pb->faceOffsets[faceB+1]; i++) {
		int idx = pb->faceIndices[i];
		pointsB[countB++] = tin_fwtrf_point(tb, pb->vertices[idx]);
	}

	// Clip the faces against each other to create a collision manifold.
	Tin_Vec3 manifold[32];
	int count = tin_clip_manifolds(pointsA, countA, pointsB, countB, refNormal, manifold);

	// Simplify the collision manifold.
	while (count > TIN_MAX_CONTACTS) {
		count = tin_reduce_manifold(manifold, count);
	}

	tin_fill_arbiter(ta, tb, faceA, faceB, refNormal, refBase, manifold, count, arbiter);

	return count;
}

Tin_Scalar
tin_dot_v12(const Tin_Scalar a[12], const Tin_Scalar b[12])
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

void
tin_apply_impulse(Tin_Constraint1D *constraint, Tin_Scalar *velocity1, Tin_Scalar *velocity2)
{
#if TIN_HAS_SSE2 // and Tin_Scalar == float ...
	/* Prepare 12D velocity vector for SIMD operations */
	__m128 vel1XY = _mm_loadu_ps(&velocity1[0]);
	__m128 vel1Z  = _mm_loadu_ps(&velocity1[4]);
	__m128 vel2XY = _mm_loadu_ps(&velocity2[0]);
	__m128 vel2Z  = _mm_loadu_ps(&velocity2[4]);
	__m128 velX   = _mm_unpacklo_ps(vel1XY, vel2XY); // V1X V2X W1X W2X
	__m128 velY   = _mm_unpackhi_ps(vel1XY, vel2XY); // V1Y V2Y W1Y W2Y
	__m128 velZ   = _mm_unpacklo_ps(vel1Z,  vel2Z);  // V1Z V2Z W1Z W2Z

	/* Apply impulse */
	__m128 invMassJacobianX = _mm_loadu_ps(&constraint->invMassJacobian[0]);
	__m128 invMassJacobianY = _mm_loadu_ps(&constraint->invMassJacobian[4]);
	__m128 invMassJacobianZ = _mm_loadu_ps(&constraint->invMassJacobian[8]);
	__m128 magnitude4 = _mm_set1_ps(constraint->accumMagnitude);
	velX = _mm_add_ps(velX, _mm_mul_ps(magnitude4, invMassJacobianX));
	velY = _mm_add_ps(velY, _mm_mul_ps(magnitude4, invMassJacobianY));
	velZ = _mm_add_ps(velZ, _mm_mul_ps(magnitude4, invMassJacobianZ));

	/* Deinterleave velocities again */
	vel1XY = _mm_shuffle_ps(velX, velY,  _MM_SHUFFLE(2, 0, 2, 0));
	vel1Z  = _mm_shuffle_ps(velZ, vel1Z, _MM_SHUFFLE(3, 2, 2, 0));
	vel2XY = _mm_shuffle_ps(velX, velY,  _MM_SHUFFLE(3, 1, 3, 1));
	vel2Z  = _mm_shuffle_ps(velZ, vel2Z, _MM_SHUFFLE(3, 2, 3, 1));

	/* Write back changed velocities */
	_mm_storeu_ps(&velocity1[0], vel1XY);
	_mm_storeu_ps(&velocity1[4], vel1Z);
	_mm_storeu_ps(&velocity2[0], vel2XY);
	_mm_storeu_ps(&velocity2[4], vel2Z);
#else
	// TODO Non-SIMD tin_apply_impulse()
#error "Non-SIMD tin_apply_impulse() is not implemented yet"
#endif
}

void
tin_enforce_constraint_1d(Tin_Constraint1D *constraint, Tin_Scalar *velocity1, Tin_Scalar *velocity2)
{
#if TIN_HAS_SSE2 // TODO and Tin_Scalar == float ...
	/* Prepare 12D velocity vector for SIMD operations */
	__m128 vel1XY = _mm_loadu_ps(&velocity1[0]);
	__m128 vel1Z  = _mm_loadu_ps(&velocity1[4]);
	__m128 vel2XY = _mm_loadu_ps(&velocity2[0]);
	__m128 vel2Z  = _mm_loadu_ps(&velocity2[4]);
	__m128 velX   = _mm_unpacklo_ps(vel1XY, vel2XY); // V1X V2X W1X W2X
	__m128 velY   = _mm_unpackhi_ps(vel1XY, vel2XY); // V1Y V2Y W1Y W2Y
	__m128 velZ   = _mm_unpacklo_ps(vel1Z,  vel2Z);  // V1Z V2Z W1Z W2Z

	/* Solve for magnitude */
	__m128 jacobianX = _mm_loadu_ps(&constraint->jacobian[0]);
	__m128 jacobianY = _mm_loadu_ps(&constraint->jacobian[4]);
	__m128 jacobianZ = _mm_loadu_ps(&constraint->jacobian[8]);
	__m128 dotX = _mm_mul_ps(jacobianX, velX);
	__m128 dotY = _mm_mul_ps(jacobianY, velY);
	__m128 dotZ = _mm_mul_ps(jacobianZ, velZ);
	__m128 dot4 = _mm_add_ps(_mm_add_ps(dotX, dotY), dotZ);
	alignas(16) float dotF[4];
	_mm_store_ps(dotF, dot4);
	float dot = (dotF[0] + dotF[1]) + (dotF[2] + dotF[3]);
	Tin_Scalar magnitude = (constraint->bias - dot) * constraint->effectiveMass;

	/* Clamp magnitude to fulfill inequality */
	Tin_Scalar accumMagnitude = constraint->accumMagnitude;
	accumMagnitude += magnitude;
	accumMagnitude = MAX(accumMagnitude, constraint->minMagnitude);
	accumMagnitude = MIN(accumMagnitude, constraint->maxMagnitude);
	magnitude = accumMagnitude - constraint->accumMagnitude;
	constraint->accumMagnitude = accumMagnitude;

	/* Apply impulse */
	__m128 invMassJacobianX = _mm_loadu_ps(&constraint->invMassJacobian[0]);
	__m128 invMassJacobianY = _mm_loadu_ps(&constraint->invMassJacobian[4]);
	__m128 invMassJacobianZ = _mm_loadu_ps(&constraint->invMassJacobian[8]);
	__m128 magnitude4 = _mm_set1_ps(magnitude);
	velX = _mm_add_ps(velX, _mm_mul_ps(magnitude4, invMassJacobianX));
	velY = _mm_add_ps(velY, _mm_mul_ps(magnitude4, invMassJacobianY));
	velZ = _mm_add_ps(velZ, _mm_mul_ps(magnitude4, invMassJacobianZ));

	/* Deinterleave velocities again */
	vel1XY = _mm_shuffle_ps(velX, velY,  _MM_SHUFFLE(2, 0, 2, 0));
	vel1Z  = _mm_shuffle_ps(velZ, vel1Z, _MM_SHUFFLE(3, 2, 2, 0));
	vel2XY = _mm_shuffle_ps(velX, velY,  _MM_SHUFFLE(3, 1, 3, 1));
	vel2Z  = _mm_shuffle_ps(velZ, vel2Z, _MM_SHUFFLE(3, 2, 3, 1));

	/* Write back changed velocities */
	_mm_storeu_ps(&velocity1[0], vel1XY);
	_mm_storeu_ps(&velocity1[4], vel1Z);
	_mm_storeu_ps(&velocity2[0], vel2XY);
	_mm_storeu_ps(&velocity2[4], vel2Z);
#else
	// TODO Non-SIMD tin_enforce_jacobian_fast()
#error "Non-SIMD tin_enforce_jacobian_fast() is not implemented yet"
#endif
}

void
tin_build_constraint_1d(Tin_Constraint1D *constraint,
	Tin_Scalar invMass1, Tin_Scalar invInertia1[9],
	Tin_Scalar invMass2, Tin_Scalar invInertia2[9],
	Tin_Vec3 posFrom1, Tin_Vec3 posFrom2,
	Tin_Vec3 axis)
{
	Tin_Vec3 r1Xaxis = tin_cross_v3(posFrom1, axis);
	Tin_Vec3 r2Xaxis = tin_cross_v3(posFrom2, axis);

	/* 12x1 Jacobian Matrix (J) */
	// SIMD friendly layout: V1X V2X W1X W2X V1Y V2Y W1Y W2Y V1Z V2Z W1Z W2Z
	constraint->jacobian[ 0] = -axis.c[0];
	constraint->jacobian[ 1] =  axis.c[0];
	constraint->jacobian[ 2] = -r1Xaxis.c[0];
	constraint->jacobian[ 3] =  r2Xaxis.c[0];
	constraint->jacobian[ 4] = -axis.c[1];
	constraint->jacobian[ 5] =  axis.c[1];
	constraint->jacobian[ 6] = -r1Xaxis.c[1];
	constraint->jacobian[ 7] =  r2Xaxis.c[1];
	constraint->jacobian[ 8] = -axis.c[2];
	constraint->jacobian[ 9] =  axis.c[2];
	constraint->jacobian[10] = -r1Xaxis.c[2];
	constraint->jacobian[11] =  r2Xaxis.c[2];

	/* Effective Mass (J M^-1 J^T) */
	Tin_Vec3 angularImpulse1 = tin_m3_times_v3(invInertia1,
			TIN_VEC3(constraint->jacobian[2], constraint->jacobian[6], constraint->jacobian[10]));
	Tin_Vec3 angularImpulse2 = tin_m3_times_v3(invInertia2,
			TIN_VEC3(constraint->jacobian[3], constraint->jacobian[7], constraint->jacobian[11]));
	constraint->invMassJacobian[ 0] = invMass1 * constraint->jacobian[0];
	constraint->invMassJacobian[ 1] = invMass2 * constraint->jacobian[1];
	constraint->invMassJacobian[ 2] = angularImpulse1.c[0];
	constraint->invMassJacobian[ 3] = angularImpulse2.c[0];
	constraint->invMassJacobian[ 4] = invMass1 * constraint->jacobian[4];
	constraint->invMassJacobian[ 5] = invMass2 * constraint->jacobian[5];
	constraint->invMassJacobian[ 6] = angularImpulse1.c[1];
	constraint->invMassJacobian[ 7] = angularImpulse2.c[1];
	constraint->invMassJacobian[ 8] = invMass1 * constraint->jacobian[8];
	constraint->invMassJacobian[ 9] = invMass2 * constraint->jacobian[9];
	constraint->invMassJacobian[10] = angularImpulse1.c[2];
	constraint->invMassJacobian[11] = angularImpulse2.c[2];
	constraint->effectiveMass = 1.0 / tin_dot_v12(constraint->jacobian, constraint->invMassJacobian);

	constraint->bias = 0.0;
	constraint->minMagnitude = -INFINITY;
	constraint->maxMagnitude =  INFINITY;
	constraint->accumMagnitude = 0.0;
}

void
tin_arbiter_prestep(Tin_Scene *scene, Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6], Tin_Scalar invDt)
{
	const Tin_Scalar allowedPenetration = 0.01;
	const Tin_Scalar biasFactor = 0.1;

	Tin_Scalar *vel1 = velocities[arbiter->bodyID1];
	Tin_Scalar *vel2 = velocities[arbiter->bodyID2];
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

	Tin_Body *body1 = &scene->bodyTable[arbiter->bodyID1];
	Tin_Body *body2 = &scene->bodyTable[arbiter->bodyID2];

	for (int i = 0; i < arbiter->numPenetrating; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		tin_build_constraint_1d(&contact->normalConstraint,
			body1->invMass, body1->invInertia,
			body2->invMass, body2->invInertia,
			contact->posFrom1, contact->posFrom2,
			arbiter->normal);
		contact->normalConstraint.minMagnitude = 0.0;
		contact->normalConstraint.bias = -biasFactor * invDt * MIN(0.0, contact->separation + allowedPenetration);

		tin_build_constraint_1d(&contact->tangentConstraint1,
			body1->invMass, body1->invInertia,
			body2->invMass, body2->invInertia,
			contact->posFrom1, contact->posFrom2,
			arbiter->frictionDir);

		tin_build_constraint_1d(&contact->tangentConstraint2,
			body1->invMass, body1->invInertia,
			body2->invMass, body2->invInertia,
			contact->posFrom1, contact->posFrom2,
			arbiter->orthoDir);
	}
}

Tin_Scalar
tin_alignof_v12(const Tin_Scalar a[12], const Tin_Scalar b[12])
{
	Tin_Scalar aSq = tin_dot_v12(a, a);
	Tin_Scalar bSq = tin_dot_v12(b, b);
	Tin_Scalar abLen = sqrt(aSq * bSq);
	Tin_Scalar proj = tin_dot_v12(a, b);
	Tin_Scalar projScaled = abLen <= TIN_EPSILON ? 0.0 : (proj / abLen);
	return MAX(projScaled, 0.0);
}

void
tin_arbiter_warm_start(Tin_Scalar (*velocities)[6], Tin_Arbiter *arbiter, const Tin_Arbiter *oldArbiter)
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
				Tin_Scalar normalProj = tin_alignof_v12(contact->normalConstraint.jacobian, oldContact->normalConstraint.jacobian);
				contact->normalConstraint.accumMagnitude = 0.8 * normalProj * oldContact->normalConstraint.accumMagnitude;
				Tin_Scalar tangentProj1 = tin_alignof_v12(contact->tangentConstraint1.jacobian, oldContact->tangentConstraint1.jacobian);
				contact->tangentConstraint1.accumMagnitude = 0.8 * tangentProj1 * oldContact->tangentConstraint1.accumMagnitude;
				Tin_Scalar tangentProj2 = tin_alignof_v12(contact->tangentConstraint2.jacobian, oldContact->tangentConstraint2.jacobian);
				contact->tangentConstraint2.accumMagnitude = 0.8 * tangentProj2 * oldContact->tangentConstraint2.accumMagnitude;

				tin_apply_impulse(&contact->normalConstraint, velocities[arbiter->bodyID1], velocities[arbiter->bodyID2]);
				tin_apply_impulse(&contact->tangentConstraint1, velocities[arbiter->bodyID1], velocities[arbiter->bodyID2]);
				tin_apply_impulse(&contact->tangentConstraint2, velocities[arbiter->bodyID1], velocities[arbiter->bodyID2]);
				break;
			}
		}
	}
}

void
tin_arbiter_apply_separation(Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6])
{
	for (int idx = 0; idx < arbiter->numPenetrating; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];
		tin_enforce_constraint_1d(&contact->normalConstraint, velocities[arbiter->bodyID1], velocities[arbiter->bodyID2]);
	}
}

void
tin_arbiter_apply_friction(Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6])
{
	const Tin_Scalar friction = 0.4;

	for (int idx = 0; idx < arbiter->numPenetrating; idx++) {
		Tin_Contact *contact = &arbiter->contacts[idx];
		Tin_Scalar frictionLimit = contact->normalConstraint.accumMagnitude * friction;
		contact->tangentConstraint1.minMagnitude = -frictionLimit;
		contact->tangentConstraint1.maxMagnitude =  frictionLimit;
		contact->tangentConstraint2.minMagnitude = -frictionLimit;
		contact->tangentConstraint2.maxMagnitude =  frictionLimit;
		tin_enforce_constraint_1d(&contact->tangentConstraint1, velocities[arbiter->bodyID1], velocities[arbiter->bodyID2]);
		tin_enforce_constraint_1d(&contact->tangentConstraint2, velocities[arbiter->bodyID1], velocities[arbiter->bodyID2]);
	}
}

/* === Island Detection & Freezing === :island: */

Tin_BodyID
tin_island_find(Tin_Scene *scene, Tin_BodyID bodyID)
{
	Tin_BodyID islandID = bodyID;
	while (scene->bodyTable[islandID].islandID != islandID) {
		islandID = scene->bodyTable[islandID].islandID;
	}
	if (scene->bodyTable[bodyID].islandID != islandID) {
		scene->bodyTable[bodyID].islandID = islandID;
	}
	return islandID;
}

void
tin_island_union(Tin_Scene *scene, Tin_BodyID bodyID1, Tin_BodyID bodyID2)
{
	bodyID1 = tin_island_find(scene, bodyID1);
	bodyID2 = tin_island_find(scene, bodyID2);
	scene->bodyTable[bodyID2].islandID = bodyID1;
}

void
tin_build_islands(Tin_Scene *scene)
{
	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];
		body->islandID = i;
		body->islandStable = true;
	}

	for (size_t i = 0; i < scene->numArbiters; i++) {
		Tin_Arbiter *arbiter = &scene->arbiters[i];
		Tin_Body *body1 = &scene->bodyTable[arbiter->bodyID1];
		Tin_Body *body2 = &scene->bodyTable[arbiter->bodyID2];
		if (body1->invMass != 0.0 && body2->invMass != 0.0) {
			tin_island_union(scene, arbiter->bodyID1, arbiter->bodyID2);
		}
	}

	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];
		if (body->restCounter >= 5 || body->invMass == 0.0) continue;
		Tin_BodyID islandID = tin_island_find(scene, i);
		scene->bodyTable[islandID].islandStable = false;
	}
}

/* === Scenes / Worlds === :scene: */

void
tin_scene_update(Tin_Scene *scene)
{
	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];

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
		tin_arbiter_prestep(scene, &scene->arbiters[i], velocities, invDt);
#if 1
		void *payload;
		if (tin_find_pair(&scene->contactCache, scene->arbiters[i].bodyID1, scene->arbiters[i].bodyID2, &payload)) {
			tin_arbiter_warm_start(velocities, &scene->arbiters[i], payload);
		}
#endif
	}
}

void
tin_scene_step(Tin_Scene *scene, Tin_Scalar (*velocities)[6])
{
	for (size_t i = 0; i < scene->numArbiters; i++) {
		tin_arbiter_apply_separation(&scene->arbiters[i], velocities);
	}
	for (size_t i = 0; i < scene->numArbiters; i++) {
		tin_arbiter_apply_friction(&scene->arbiters[i], velocities);
	}
}

void
tin_integrate(Tin_Scene *scene, Tin_Scalar dt)
{
	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];

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
tin_narrowphase(Tin_Scene *scene)
{
	for (size_t i = 0; i < scene->numArbiters; i++) {
		Tin_Arbiter *arbiter = &scene->arbiters[i];
		Tin_Body *body1 = &scene->bodyTable[arbiter->bodyID1];
		Tin_Body *body2 = &scene->bodyTable[arbiter->bodyID2];
		arbiter->numContacts = tin_polytope_collide(
				&body1->shape->polytope, &body1->transform,
				&body2->shape->polytope, &body2->transform, &scene->arbiters[i]);
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
	tin_sweep_prune_update(scene->sweepPrune);
	Tin_BodyID *broad;
	size_t numBroad;
	tin_sweep_prune(scene->sweepPrune, &broad, &numBroad);
	scene->numArbiters = 0;
	for (size_t i = 0; i < numBroad; i += 2) {
		if (scene->numArbiters == scene->capArbiters) {
			scene->capArbiters = scene->capArbiters ? 2*scene->capArbiters : 16;
			scene->arbiters = realloc(scene->arbiters, scene->capArbiters * sizeof *scene->arbiters);
		}
		scene->arbiters[scene->numArbiters++] = (Tin_Arbiter){
			.bodyID1 = broad[i],
			.bodyID2 = broad[i+1],
		};
	}
	free(broad);

	tin_build_islands(scene);
	/* Apply Gravity */
	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];
		if (!scene->bodyTable[tin_island_find(scene, i)].islandStable) {
			body->velocity = tin_saxpy_v3(dt, TIN_VEC3(0.0, -9.0, 0.0), body->velocity);
		}
	}

	/* Filter out collisions within stable islands */
	size_t oldNumArbiters = scene->numArbiters;
	scene->numArbiters = 0;
	for (size_t i = 0; i < oldNumArbiters; i++) {
		Tin_Arbiter *arbiter = &scene->arbiters[i];
		if (scene->bodyTable[tin_island_find(scene, arbiter->bodyID1)].islandStable) {
			if (scene->bodyTable[tin_island_find(scene, arbiter->bodyID2)].islandStable) {
				continue;
			}
		}
		scene->arbiters[scene->numArbiters++] = *arbiter;
	}
	printf("#collisions = %zu\n", scene->numArbiters);

	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[1] += stopTime - startTime;

	startTime = gettime ? gettime() : 0.0;
	tin_narrowphase(scene);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[2] += stopTime - startTime;
	
	// We keep some padding at the end so we can more easily load values into SIMD registers.
	Tin_Scalar (*velocities)[6] = calloc((scene->bodyTableCapac + 1), sizeof (Tin_Scalar[6]));
	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];

		velocities[i][0] = body->velocity.c[0];
		velocities[i][1] = body->angularVelocity.c[0];
		velocities[i][2] = body->velocity.c[1];
		velocities[i][3] = body->angularVelocity.c[1];
		velocities[i][4] = body->velocity.c[2];
		velocities[i][5] = body->angularVelocity.c[2];
	}
	
	startTime = gettime ? gettime() : 0.0;
	tin_scene_prestep(scene, velocities, invDt);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[3] += stopTime - startTime;

	startTime = gettime ? gettime() : 0.0;
	for (int iter = 0; iter < 16; iter++) {
		tin_scene_step(scene, velocities);
	}
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[4] += stopTime - startTime;

	for (Tin_BodyID i = 0; i < scene->bodyTableCapac; i++) {
		if (!scene->bodyOccupied[i]) continue;
		Tin_Body *body = &scene->bodyTable[i];

		body->velocity.c[0] = velocities[i][0];
		body->angularVelocity.c[0] = velocities[i][1];
		body->velocity.c[1] = velocities[i][2];
		body->angularVelocity.c[1] = velocities[i][3];
		body->velocity.c[2] = velocities[i][4];
		body->angularVelocity.c[2] = velocities[i][5];
	}
	free(velocities);

	startTime = gettime ? gettime() : 0.0;
	// Cache this frames contacts for the next frame
	tin_reset_pairtable(&scene->contactCache);
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
			arbiter->bodyID1, arbiter->bodyID2,
			arbiter);
	}
	tin_integrate(scene, dt);
	stopTime = gettime ? gettime() : 0.0;
	if (timings) timings[5] += stopTime - startTime;
}

Tin_BodyID
tin_add_body(Tin_Scene *scene, const Tin_Shape *shape, Tin_Scalar invMass)
{
	if (scene->freeBodyID == scene->bodyTableCapac) {
		scene->bodyTableCapac *= 2;
		if (!scene->bodyTableCapac) scene->bodyTableCapac = 16;
		scene->bodyTable = realloc(scene->bodyTable,
			scene->bodyTableCapac * sizeof *scene->bodyTable);
		scene->bodyOccupied = realloc(scene->bodyOccupied,
			scene->bodyTableCapac * sizeof *scene->bodyOccupied);
		for (Tin_BodyID i = scene->freeBodyID; i < scene->bodyTableCapac; i++) {
			scene->bodyOccupied[i] = false;
			Tin_BodyID *nextField = (Tin_BodyID *)&scene->bodyTable[i];
			*nextField = i + 1;
		}
	}

	Tin_BodyID bodyID = scene->freeBodyID;
	scene->bodyOccupied[bodyID] = true;
	Tin_BodyID *nextField = (Tin_BodyID *)&scene->bodyTable[bodyID];
	scene->freeBodyID = *nextField;

	Tin_Body *body = &scene->bodyTable[bodyID];
	memset(body, 0, sizeof *body);
	body->transform = (Tin_Transform) {
		{ 1, 0, 0, 0, 1, 0, 0, 0, 1 },
		TIN_VEC3(0, 0, 0),
		1.0,
	};
	body->shape = shape;
	body->invMass = invMass;
	tin_sweep_prune_add_body(scene->sweepPrune, bodyID);
	return bodyID;
}

void
tin_delete_body(Tin_Scene *scene, Tin_BodyID bodyID)
{
	if (bodyID >= scene->bodyTableCapac) return;
	if (!scene->bodyOccupied[bodyID]) return;
	scene->bodyOccupied[bodyID] = false;
	tin_sweep_prune_delete_body(scene->sweepPrune, bodyID);
	Tin_BodyID *nextField = (Tin_BodyID *)&scene->bodyTable[bodyID];
	*nextField = scene->freeBodyID;
	scene->freeBodyID = bodyID;
}
