/* tincan physics engine
 *
 * Copyright 2021-2025 Thomas Oltmann
 * All rights reserved
 */

#ifndef TINCAN_PHYSICS_H
#define TINCAN_PHYSICS_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h> /* for offsetof */
#include <float.h>

/* === Linked List === :list: */

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

/* === Customizable Allocators === :alloc: */

typedef struct {
	void *userPointer;
	void *(*alloc)(void *userPointer);
	void (*free)(void *userPointer, void *memoryPointer);
} Tin_Allocator;

/* === Vector Math === :vec: */

#define TIN_EPSILON FLT_EPSILON
typedef float Tin_Scalar;

struct Tin_Vec3 { Tin_Scalar c[3]; };
typedef struct Tin_Vec3 Tin_Vec3;

#define TIN_VEC3(x, y, z) (Tin_Vec3){{x,y,z}}

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

void tin_axis_angle_to_matrix(Tin_Vec3 axis, Tin_Scalar angle, Tin_Scalar matrix[3*3]);
void tin_m3_times_m3(Tin_Scalar result[3*3], const Tin_Scalar matrixA[3*3], const Tin_Scalar matrixB[3*3]);
void tin_m3_times_m3_transposed(Tin_Scalar result[3*3], const Tin_Scalar matrixA[3*3], const Tin_Scalar matrixB[3*3]);
Tin_Vec3 tin_m3_times_v3(const Tin_Scalar matrix[3*3], Tin_Vec3 vector);

Tin_Vec3   tin_gram_schmidt(Tin_Vec3 fixed, Tin_Vec3 var);
Tin_Scalar tin_prlgram_area(Tin_Vec3 e1, Tin_Vec3 e2);

/* === Transforms === :trf: */

typedef struct {
	Tin_Scalar rotation[3*3];
	Tin_Vec3   translation;
	Tin_Scalar scale;
} Tin_Transform;

Tin_Vec3 tin_fwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec);
Tin_Vec3 tin_fwtrf_dir  (const Tin_Transform *transform, Tin_Vec3 vec);
Tin_Vec3 tin_bwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec);
Tin_Vec3 tin_bwtrf_dir  (const Tin_Transform *transform, Tin_Vec3 vec);

/* === Polytopes === :poly: */

typedef struct {
	Tin_Vec3  *vertices;
	int       *faceIndices;
	int       *faceOffsets;
	Tin_Vec3  *faceNormals;
	int        numVertices;
	int        numFaces;
} Tin_Polytope;

Tin_Vec3 tin_polytope_support(const void *geometry, Tin_Vec3 dir);

/* === Shapes === :shape: */

#define TIN_SPHERE   's'
#define TIN_POLYTOPE 'p'

typedef struct {
	Tin_Vec3   invInertia;
	Tin_Scalar radius;
	int        kind;
	union {
		Tin_Polytope polytope;
	};
} Tin_Shape;

void tin_shape_aabb(const Tin_Shape *shape, const Tin_Transform *transform, Tin_Vec3 *aabbMin, Tin_Vec3 *aabbMax);

/* === Rigid Bodies === :body: */

typedef struct Tin_Body Tin_Body;
struct Tin_Body {
	Tin_List      node;
	Tin_Transform transform;
	const Tin_Shape *shape;
	Tin_Scalar    invMass;
	Tin_Vec3      velocity;
	Tin_Vec3      angularVelocity;
	Tin_Scalar    invInertia[3*3];
	Tin_Vec3      aabbMin;
	Tin_Vec3      aabbMax;
	Tin_Body     *island;
	int           restCounter;
	bool          islandStable;
};

/* === Minkowski Sum === :mink: */
/* (difference of transformed convex polytopes) */
typedef struct {
	const Tin_Polytope  *polytope1;
	const Tin_Transform *transform1;
	const Tin_Polytope  *polytope2;
	const Tin_Transform *transform2;
} Tin_Polysum;

Tin_Vec3 tin_polysum_support(const void *geometry, Tin_Vec3 dir);

/* === Minkowski Portal Refinement === :mpr: */

typedef struct {
	Tin_Vec3 origin;
	Tin_Vec3 dir;
} Tin_Ray;

typedef struct {
	Tin_Vec3 a, b, c;
	Tin_Vec3 normal;
} Tin_Portal;

typedef Tin_Vec3 Tin_SupportFunc(const void *, Tin_Vec3);

int tin_construct_portal(const void *geometry, Tin_SupportFunc support,
	const Tin_Ray *r, Tin_Portal *p);
void tin_refine_portal  (const void *geometry, Tin_SupportFunc support,
	const Tin_Ray *r, Tin_Portal *p);

/* === Contact Points === :contact: */

#define TIN_MAX_CONTACTS 4

typedef struct {
	Tin_Scalar jacobian[12];
	Tin_Vec3   angularImpulse1;
	Tin_Vec3   angularImpulse2;
	Tin_Scalar effectiveMass;
	Tin_Scalar bias;
	Tin_Scalar minMagnitude;
	Tin_Scalar maxMagnitude;
	Tin_Scalar accumMagnitude;
} Tin_Constraint1D;

typedef struct {
	/* Vector from origin of body to the contact point on the body, in world-space. */
	Tin_Vec3   posFrom1;
	Tin_Vec3   posFrom2;

	Tin_Scalar separation;

	Tin_Constraint1D normalConstraint;
	Tin_Constraint1D tangentConstraint1;
	Tin_Constraint1D tangentConstraint2;
} Tin_Contact;

/*
typedef struct {
	Tin_Vec3 posFrom1;
	Tin_Vec3 posFrom2;
	Tin_Scalar magnitudeAccums[3];
} Tin_CachedContactPoint;

typedef struct {
	int face1;
	int face2;
	Tin_CachedContactPoint points[TIN_MAX_CONTACTS];
	int numPoints;
} Tin_CachedContact;
*/

typedef struct {
	Tin_Body *body1;
	Tin_Body *body2;
	int body1Idx;
	int body2Idx;
	Tin_Scalar body1InvMass;
	Tin_Scalar body2InvMass;
	int face1;
	int face2;
	Tin_Vec3 normal;
	Tin_Vec3 frictionDir;
	Tin_Vec3 orthoDir;
	int numContacts;
	int numPenetrating; // the penetrating contacts (separation < 0) come first
	Tin_Contact contacts[TIN_MAX_CONTACTS];
} Tin_Arbiter;

void tin_arbiter_prestep(Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6], Tin_Scalar invDt);

/* === Pair-Indexed Hashtable === :pair: */

typedef struct {
	size_t   elemLow;
	size_t   elemHigh;
	void    *payload;
	uint32_t hash;
	bool     occupied;
} Tin_PairTableSlot;

typedef struct {
	Tin_PairTableSlot *slots;
	size_t count;
	size_t capac;
} Tin_PairTable;

void tin_create_pairtable(Tin_PairTable *table);
void tin_destroy_pairtable(Tin_PairTable *table);
void tin_reset_pairtable(Tin_PairTable *table);
void tin_resize_pairtable(Tin_PairTable *table, size_t newCapac);
void tin_order_pair(size_t *elemLow, size_t *elemHigh);
uint32_t tin_hash_pair(size_t elemLow, size_t elemHigh);
size_t tin_fold_hash(size_t capac, uint32_t hash);
size_t tin_pairtable_index(Tin_PairTable *table, size_t elemLow, size_t elemHigh);
bool tin_find_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void **payloadOut);
void tin_insert_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void *payload);
void tin_delete_pair(Tin_PairTable *table, size_t elemA, size_t elemB);

/* === Island Detection & Freezing === :island: */

typedef struct Tin_Scene Tin_Scene;

Tin_Body *tin_island_find(Tin_Body *body);
void tin_island_union(Tin_Body *body1, Tin_Body *body2);
void tin_build_islands(Tin_Scene *scene);

/* === Scenes / Worlds === :scene: */

typedef struct Tin_Scene {
	Tin_List bodies;
	Tin_Allocator bodyAllocator;
	Tin_PairTable contactCache;
	Tin_Arbiter *arbiters;
	size_t numArbiters;
	size_t capArbiters;
	Tin_Arbiter *oldArbiters;
	size_t numOldArbiters;
	size_t capOldArbiters;
} Tin_Scene;

void tin_scene_update(Tin_Scene *scene);
void tin_scene_prestep(Tin_Scene *scene, Tin_Scalar (*velocities)[6], Tin_Scalar invDt);
void tin_scene_step(Tin_Scene *scene, Tin_Scalar (*velocities)[6]);
void tin_integrate(Tin_Scene *scene, Tin_Scalar dt);
void tin_broadphase(Tin_Scene *scene);
void tin_simulate(Tin_Scene *scene, Tin_Scalar dt, double (*gettime)(), double timings[6]);

Tin_Body *tin_add_body(Tin_Scene *scene, const Tin_Shape *shape, Tin_Scalar invMass);

#endif
