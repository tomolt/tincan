#include "tincan.h"

#include <string.h>
#include <tgmath.h>
#include <assert.h>

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

