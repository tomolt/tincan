#ifndef TINCAN_PHYSICS_H
#define TINCAN_PHYSICS_H

#include <stdint.h>

typedef float tin_scalar;

/* 3D Vectors */

struct tin_vec3 { tin_scalar c[3]; };
typedef struct tin_vec3 tin_vec3;

tin_vec3 tin_neg_v3  (tin_vec3 x);
tin_vec3 tin_add_v3  (tin_vec3 a, tin_vec3 b);
tin_vec3 tin_sub_v3  (tin_vec3 a, tin_vec3 b);
tin_vec3 tin_scale_v3(tin_scalar a, tin_vec3 x);
tin_vec3 tin_saxpy_v3(tin_scalar a, tin_vec3 x, tin_vec3 y);
tin_vec3 tin_cross_v3(tin_vec3 a, tin_vec3 b);
tin_scalar tin_dot_v3(tin_vec3 a, tin_vec3 b);
tin_vec3 tin_normalize_v3(tin_vec3 v);

/* Quaternions */

struct tin_quat { tin_vec3 v; tin_scalar s; };
typedef struct tin_quat tin_quat;

tin_quat tin_make_qt (tin_vec3 axis, tin_scalar angle);
tin_quat tin_mul_qt  (tin_quat a, tin_quat b);
tin_vec3 tin_apply_qt(tin_quat q, tin_vec3 v);
tin_quat tin_conjugate_qt(tin_quat q);

/* Geometric Functions */

tin_vec3   tin_gram_schmidt(tin_vec3 fixed, tin_vec3 var);
tin_scalar tin_prlgram_area(tin_vec3 e1, tin_vec3 e2);

/* Transforms */

typedef struct {
	tin_quat   rotation;
	tin_vec3   translation;
	tin_scalar scale;
} tin_transform;

/* Shapes */

#define TIN_SPHERE 's'
#define TIN_CONVEX 'c'

/* Rigid Bodies */

typedef struct {
	tin_transform transform;
	const void   *shape_params;
	int           shape;
} tin_body;

/* Polytopes */

typedef struct {
	tin_vec3 *vertices;
	int       num_vertices;
} tin_polytope;

tin_vec3 tin_polytope_support(const tin_polytope *p, tin_vec3 dir);

/* Minkowski sum (difference) of transformed convex polytopes */
typedef struct {
	const tin_polytope  *former_polytope;
	const tin_transform *former_transform;
	const tin_polytope  *latter_polytope;
	const tin_transform *latter_transform;
} tin_polysum;

/* A point in a polytope sum */
typedef struct {
	tin_vec3 abs;
	tin_vec3 rel_former;
	tin_vec3 rel_latter;
} tin_pspoint;

void tin_polysum_support(const tin_polysum *s, tin_vec3 dir, tin_pspoint *sup);

/* MPR */

typedef struct {
	tin_vec3 origin;
	tin_vec3 dir;
} tin_ray;

typedef struct {
	tin_pspoint a, b, c;
	tin_vec3    normal;
} tin_portal;

int  tin_construct_portal(const tin_polysum *ps, const tin_ray *r, tin_portal *p);
void tin_refine_portal   (const tin_polysum *ps, const tin_ray *r, tin_portal *p);

#if 0

typedef struct {
	tin_vec3   former_rel;
	tin_vec3   latter_rel;
	tin_vec3   normal;
	tin_scalar depth;
} tin_contact;

void minkowski_support(const Body *b[2], tin_vec3 dir, SharedPoint *sup);
int  mpr_intersect    (const Body *b[2], const Ray *ray, Contact *contact);

#endif

#endif
