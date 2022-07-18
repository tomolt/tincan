#ifndef TINCAN_PHYSICS_H
#define TINCAN_PHYSICS_H

#include <stdint.h>
#include <float.h>

#define TIN_EPSILON FLT_EPSILON
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
tin_vec3 tin_hadamard_v3(tin_vec3 a, tin_vec3 b);

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

tin_vec3 tin_fwtrf_point(const tin_transform *transform, tin_vec3 vec);
tin_vec3 tin_fwtrf_dir  (const tin_transform *transform, tin_vec3 vec);
tin_vec3 tin_bwtrf_point(const tin_transform *transform, tin_vec3 vec);
tin_vec3 tin_bwtrf_dir  (const tin_transform *transform, tin_vec3 vec);

/* Shapes */

#define TIN_SPHERE 's'
#define TIN_CONVEX 'c'

/* Rigid Bodies */

typedef struct {
	tin_transform transform;
	const void   *shape_params;
	int           shape;
	tin_scalar    inv_mass;
	tin_vec3      inv_inertia;
	tin_vec3      velocity;
	tin_vec3      angular_velocity;
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

/* Contact Points */

typedef struct {
	tin_vec3   rel1;
	tin_vec3   rel2;
	tin_vec3   normal;
	tin_scalar base_stretch;

	tin_vec3   position;
	tin_scalar separation;
	
	tin_scalar normal_mass;
	tin_scalar bias;
} tin_contact;

int tin_polytope_collide(
	const tin_polytope *pa, const tin_transform *ta,
	const tin_polytope *pb, const tin_transform *tb,
	tin_contact *contact);

#define TIN_MAX_CONTACTS 5

typedef struct {
	tin_body *body1;
	tin_body *body2;
	int num_contacts;
	tin_contact contacts[TIN_MAX_CONTACTS];
} tin_arbiter;

void tin_arbiter_update(tin_arbiter *a);
void tin_arbiter_add_contact(tin_arbiter *a, tin_contact contact);
void tin_arbiter_prestep(tin_arbiter *a, tin_scalar inv_dt);
void tin_arbiter_apply_impulse(tin_arbiter *a, tin_scalar inv_dt);

#endif
