#ifndef TINCAN_PHYSICS_H
#define TINCAN_PHYSICS_H

#include <stdint.h>
#include <stddef.h> /* for offsetof */
#include <float.h>

/* Linked List */

typedef struct Tin_List Tin_List;

struct Tin_List {
	Tin_List *next;
	Tin_List *prev;
};

#define TIN_LIST_LINK(node1, node2) do { (node1).next = &(node2); (node2).prev = &(node1); } while (0)
#define TIN_LIST_INIT(list) TIN_LIST_LINK(list, list)
#define TIN_CONTAINER_OF(ptr, type, member) ((type *)((char *)(ptr) - offsetof(type, member)))
#define TIN_FOR_RANGE(elem, start, end, type, member)						\
	type *elem;									\
	Tin_List *_node;								\
	for (elem = TIN_CONTAINER_OF(_node = (start).next, type, member);		\
	     _node != &(end);								\
	     elem = TIN_CONTAINER_OF(_node = _node->next, type, member))
#define TIN_FOR_EACH(elem, list, type, member) TIN_FOR_RANGE(elem, list, list, type, member)

/* Customizable Allocators */

typedef struct {
	void *userPointer;
	void *(*alloc)(void *userPointer);
	void (*free)(void *userPointer, void *memoryPointer);
} Tin_Allocator;

/* 3D Vectors */

#define TIN_EPSILON FLT_EPSILON
typedef float Tin_Scalar;

struct Tin_Vec3 { Tin_Scalar c[3]; };
typedef struct Tin_Vec3 Tin_Vec3;

Tin_Vec3 tin_neg_v3  (Tin_Vec3 x);
Tin_Vec3 tin_add_v3  (Tin_Vec3 a, Tin_Vec3 b);
Tin_Vec3 tin_sub_v3  (Tin_Vec3 a, Tin_Vec3 b);
Tin_Vec3 tin_scale_v3(Tin_Scalar a, Tin_Vec3 x);
Tin_Vec3 tin_saxpy_v3(Tin_Scalar a, Tin_Vec3 x, Tin_Vec3 y);
Tin_Vec3 tin_cross_v3(Tin_Vec3 a, Tin_Vec3 b);
Tin_Scalar tin_dot_v3(Tin_Vec3 a, Tin_Vec3 b);
Tin_Scalar tin_length_v3(Tin_Vec3 v);
Tin_Vec3 tin_normalize_v3(Tin_Vec3 v);
Tin_Vec3 tin_hadamard_v3(Tin_Vec3 a, Tin_Vec3 b);

/* Quaternions */

struct Tin_Quat { Tin_Vec3 v; Tin_Scalar s; };
typedef struct Tin_Quat Tin_Quat;

Tin_Quat tin_make_qt (Tin_Vec3 axis, Tin_Scalar angle);
Tin_Quat tin_mul_qt  (Tin_Quat a, Tin_Quat b);
Tin_Vec3 tin_apply_qt(Tin_Quat q, Tin_Vec3 v);
Tin_Quat tin_conjugate_qt(Tin_Quat q);

/* Geometric Functions */

Tin_Vec3   tin_gram_schmidt(Tin_Vec3 fixed, Tin_Vec3 var);
Tin_Scalar tin_prlgram_area(Tin_Vec3 e1, Tin_Vec3 e2);

/* Transforms */

typedef struct {
	Tin_Quat   rotation;
	Tin_Vec3   translation;
	Tin_Scalar scale;
} Tin_Transform;

Tin_Vec3 tin_fwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec);
Tin_Vec3 tin_fwtrf_dir  (const Tin_Transform *transform, Tin_Vec3 vec);
Tin_Vec3 tin_bwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec);
Tin_Vec3 tin_bwtrf_dir  (const Tin_Transform *transform, Tin_Vec3 vec);

/* Polytopes */

typedef struct {
	Tin_Vec3  *vertices;
	int        numVertices;
} Tin_Polytope;

Tin_Vec3 tin_polytope_support(const Tin_Polytope *polytope, Tin_Vec3 dir);

/* Shapes */

#define TIN_SPHERE   's'
#define TIN_POLYTOPE 'p'

typedef struct {
	Tin_Scalar radius;
	int        kind;
	union {
		Tin_Polytope polytope;
	};
} Tin_Shape;

/* Rigid Bodies */

typedef struct {
	Tin_List      node;
	Tin_Transform transform;
	const Tin_Shape *shape;
	Tin_Scalar    invMass;
	Tin_Vec3      invInertia;
	Tin_Vec3      velocity;
	Tin_Vec3      angularVelocity;
} Tin_Body;

/* Minkowski sum (difference) of transformed convex polytopes */
typedef struct {
	const Tin_Polytope  *polytope1;
	const Tin_Transform *transform1;
	const Tin_Polytope  *polytope2;
	const Tin_Transform *transform2;
} Tin_Polysum;

/* A point in a polytope sum */
typedef struct {
	Tin_Vec3 abs;
	Tin_Vec3 relTo1;
	Tin_Vec3 relTo2;
} Tin_Pspoint;

void tin_polysum_support(const Tin_Polysum *s, Tin_Vec3 dir, Tin_Pspoint *sup);

/* MPR */

typedef struct {
	Tin_Vec3 origin;
	Tin_Vec3 dir;
} Tin_Ray;

typedef struct {
	Tin_Pspoint a, b, c;
	Tin_Vec3    normal;
} Tin_Portal;

int  tin_construct_portal(const Tin_Polysum *ps, const Tin_Ray *r, Tin_Portal *p);
void tin_refine_portal   (const Tin_Polysum *ps, const Tin_Ray *r, Tin_Portal *p);

/* Contact Points */

typedef struct {
	Tin_Vec3   rel1;
	Tin_Vec3   rel2;
	Tin_Vec3   normal;
	Tin_Scalar baseStretch;

	Tin_Vec3   position;
	Tin_Scalar separation;
	
	Tin_Scalar normalMass;
	Tin_Scalar bias;
} Tin_Contact;

int tin_polytope_collide(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Contact *contact);

#define TIN_MAX_CONTACTS 5

typedef struct {
	Tin_List node;
	Tin_Body *body1;
	Tin_Body *body2;
	int numContacts;
	Tin_Contact contacts[TIN_MAX_CONTACTS];
} Tin_Arbiter;

Tin_Vec3 tin_solve_inertia(const Tin_Body *body, Tin_Vec3 vec);
void tin_arbiter_update(Tin_Arbiter *arbiter);
void tin_arbiter_add_contact(Tin_Arbiter *arbiter, Tin_Contact contact);
void tin_arbiter_prestep(Tin_Arbiter *arbiter, Tin_Scalar invDt);
void tin_arbiter_apply_impulse(Tin_Arbiter *arbiter, Tin_Scalar invDt);

typedef struct {
	Tin_List bodies;
	Tin_List arbiters;
	Tin_Allocator bodyAllocator;
	Tin_Allocator arbiterAllocator;
} Tin_Scene;

Tin_Arbiter *tin_find_arbiter(Tin_Scene *scene, Tin_Body *body1, Tin_Body *body2);
void tin_scene_update(Tin_Scene *scene);
void tin_check_collision(Tin_Scene *scene, Tin_Body *body1, Tin_Body *body2);
void tin_scene_prestep(Tin_Scene *scene, Tin_Scalar invDt);
void tin_scene_step(Tin_Scene *scene, Tin_Scalar invDt);
void tin_integrate(Tin_Scene *scene, Tin_Scalar dt);
void tin_broadphase(Tin_Scene *scene, void (*func)(Tin_Scene *, Tin_Body *, Tin_Body *));
void tin_simulate(Tin_Scene *scene, Tin_Scalar dt);

Tin_Body *tin_add_body(Tin_Scene *scene);
Tin_Arbiter *tin_add_arbiter(Tin_Scene *scene);

#endif
