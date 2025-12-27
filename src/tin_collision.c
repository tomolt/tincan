#include "tincan.h"

#include <tgmath.h>

bool
tin_intersect_sphere_sphere(
	const Tin_Sphere *sphere1, const Tin_Transform *transform1,
	const Tin_Sphere *sphere2, const Tin_Transform *transform2,
	Tin_Arbiter *arbiter)
{
	Tin_Scalar radius1 = transform1->scale * sphere1->radius;
	Tin_Scalar radius2 = transform2->scale * sphere2->radius;
	Tin_Scalar threshold = radius1 + radius2;
	if (!(threshold > 0.0)) {
		return false;
	}
	Tin_Vec3   diff = tin_sub_v3(transform2->translation, transform1->translation);
	Tin_Scalar distSq = tin_dot_v3(diff, diff);
	bool intersect = distSq <= threshold * threshold;

	if (!intersect) {
		return false;
	}
	Tin_Scalar dist = sqrt(distSq);
	Tin_Vec3 dir;
	if (dist > 0.0) {
		dir = tin_scale_v3(1.0 / dist, diff);
	} else {
		dir = TIN_VEC3(0.0, 1.0, 0.0);
	}
	Tin_Vec3 contact = tin_saxpy_v3(radius1 / threshold, diff, transform1->translation);
	Tin_Scalar base = tin_dot_v3(dir, contact);
	tin_fill_arbiter(transform1, transform2, 0, 0, dir, base, &contact, 1, arbiter);
	return true;
}

void
tin_shape_sphere_get_aabb(const void *shape, const Tin_Transform *transform,
	Tin_Vec3 *aabbMin, Tin_Vec3 *aabbMax)
{
	const Tin_Sphere *sphere = shape;
	Tin_Vec3 center = transform->translation;
	Tin_Scalar radius = transform->scale * sphere->radius;
	Tin_Vec3 radiusVec = TIN_VEC3(radius, radius, radius);
	*aabbMin = tin_sub_v3(center, radiusVec);
	*aabbMax = tin_add_v3(center, radiusVec);
}

bool
tin_shape_sphere_intersect(const void *shape, const Tin_Transform *transform,
	const void *otherShape, const Tin_Transform *otherTransform,
	Tin_Arbiter *arbiter)
{
	const Tin_ShapeClass *const *otherClassPtr = otherShape;
	switch ((*otherClassPtr)->code) {
	case TIN_SHAPE_CODE_SPHERE:
		return tin_intersect_sphere_sphere(shape, transform, otherShape, otherTransform, arbiter);

	default:
		return false;
	}
}

const Tin_ShapeClass tin_shape_sphere = {
	.code      = TIN_SHAPE_CODE_SPHERE,
	.get_aabb  = tin_shape_sphere_get_aabb,
	.intersect = tin_shape_sphere_intersect,
};

Tin_Vec3
tin_shape_polytope_get_inv_inertia(const void *shape)
{
	const Tin_Polytope *polytope = shape;
	return polytope->invInertia;
}

void
tin_shape_polytope_get_aabb(const void *shape, const Tin_Transform *transform,
	Tin_Vec3 *aabbMin, Tin_Vec3 *aabbMax)
{
	const Tin_Polytope *polytope = shape;
	Tin_Vec3 center = transform->translation;
	Tin_Scalar radius = transform->scale * polytope->boundRadius;
	Tin_Vec3 radiusVec = TIN_VEC3(radius, radius, radius);
	*aabbMin = tin_sub_v3(center, radiusVec);
	*aabbMax = tin_add_v3(center, radiusVec);
}

bool
tin_intersect_polytope_polytope(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Arbiter *arbiter)
{
	Tin_Vec3 normal;
	if (!tin_intersect(pa, ta, pb, tb, &normal)) {
		return false;
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

	return true;
}

bool
tin_shape_polytope_intersect(const void *shape, const Tin_Transform *transform,
	const void *otherShape, const Tin_Transform *otherTransform,
	Tin_Arbiter *arbiter)
{
	const Tin_ShapeClass *const *otherClassPtr = otherShape;
	switch ((*otherClassPtr)->code) {
	case TIN_SHAPE_CODE_POLYTOPE:
		return tin_intersect_polytope_polytope(shape, transform, otherShape, otherTransform, arbiter);

	default:
		return false;
	}
}

const Tin_ShapeClass tin_shape_polytope = {
	.code            = TIN_SHAPE_CODE_POLYTOPE,
	.get_inv_inertia = tin_shape_polytope_get_inv_inertia,
	.get_aabb        = tin_shape_polytope_get_aabb,
	.intersect       = tin_shape_polytope_intersect,
};
