/* tincan physics engine
 *
 * Copyright 2021-2025 Thomas Oltmann
 * All rights reserved
 */

#ifndef TINCAN_PHYSICS_H
#define TINCAN_PHYSICS_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <float.h>

/// \defgroup allocator Customizable Allocators
/// @{

/// Customizable allocator.
typedef struct {
	void *userPointer;
	void *(*alloc)(void *userPointer);
	void (*free)(void *userPointer, void *memoryPointer);
} Tin_Allocator;

/// @}

/// \defgroup vector Vector Math
/// @{

/// An epsilon value for the Tin_Scalar type.
#define TIN_EPSILON FLT_EPSILON

/// The floating point type used to represent real-valued scalars.
typedef float Tin_Scalar;

/// Represents a real-valued (floating point) 3D vector.
typedef struct Tin_Vec3 {
	/// The cardinal axes components of the vector.
	Tin_Scalar c[3];
} Tin_Vec3;

/// Initialize a Tin_Vec3 with the given component values.
#define TIN_VEC3(x, y, z) (Tin_Vec3){{x,y,z}}

/// Invert the direction of a 3D vector.
Tin_Vec3 tin_neg_v3  (Tin_Vec3 x);
/// Add two 3D vectors together.
Tin_Vec3 tin_add_v3  (Tin_Vec3 a, Tin_Vec3 b);
/// Subtract a 3D vector from another.
Tin_Vec3 tin_sub_v3  (Tin_Vec3 a, Tin_Vec3 b);
/// Scale a 3D vector by a scalar value.
Tin_Vec3 tin_scale_v3(Tin_Scalar a, Tin_Vec3 x);
Tin_Vec3 tin_saxpy_v3(Tin_Scalar a, Tin_Vec3 x, Tin_Vec3 y);
/// Calculate the cross product of two 3D vectors.
Tin_Vec3 tin_cross_v3(Tin_Vec3 a, Tin_Vec3 b);
/// Calculate the dot product of two 3D vectors.
Tin_Scalar tin_dot_v3(Tin_Vec3 a, Tin_Vec3 b);
Tin_Scalar tin_length_v3(Tin_Vec3 v);
Tin_Vec3 tin_normalize_v3(Tin_Vec3 v);
Tin_Vec3 tin_hadamard_v3(Tin_Vec3 a, Tin_Vec3 b);

void tin_axis_angle_to_matrix(Tin_Vec3 axis, Tin_Scalar angle, Tin_Scalar matrix[3*3]);
void tin_m3_times_m3(Tin_Scalar result[3*3], const Tin_Scalar matrixA[3*3], const Tin_Scalar matrixB[3*3]);
Tin_Vec3 tin_v3_times_m3(Tin_Vec3 vector, const Tin_Scalar matrix[3*3]);
void tin_m3_times_m3_transposed(Tin_Scalar result[3*3], const Tin_Scalar matrixA[3*3], const Tin_Scalar matrixB[3*3]);
Tin_Vec3 tin_m3_times_v3(const Tin_Scalar matrix[3*3], Tin_Vec3 vector);

Tin_Vec3   tin_gram_schmidt(Tin_Vec3 fixed, Tin_Vec3 var);
Tin_Scalar tin_prlgram_area(Tin_Vec3 e1, Tin_Vec3 e2);

/// Position, orientation, and scale in 3D space.
typedef struct {
	Tin_Scalar rotation[3*3]; ///< Rotation matrix, in column-major order
	Tin_Vec3   translation; ///< Translation relative to origin
	Tin_Scalar scale; ///< Uniform scale factor
} Tin_Transform;

/// Transform a point in space.
Tin_Vec3 tin_fwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec);
/// Transform a direction vector.
Tin_Vec3 tin_fwtrf_dir  (const Tin_Transform *transform, Tin_Vec3 vec);
/// Apply inverse transform to a point in space.
Tin_Vec3 tin_bwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec);
/// Apply inverse transform to a direction vector.
Tin_Vec3 tin_bwtrf_dir  (const Tin_Transform *transform, Tin_Vec3 vec);

/// @}

/// \defgroup polytope Polytopes
/// @{

typedef struct Tin_ShapeClass Tin_ShapeClass;

typedef struct {
	const Tin_ShapeClass *vtable;
	Tin_Vec3   invInertia;
	Tin_Scalar boundRadius;
	Tin_Vec3  *vertices;
	int       *faceIndices;
	int       *faceOffsets;
	Tin_Vec3  *faceNormals;
	int        numVertices;
	int        numFaces;
} Tin_Polytope;

Tin_Vec3 tin_polytope_support(const void *geometry, Tin_Vec3 dir);

/// @}

/// \defgroup shape Collision Shapes
/// @{

#define TIN_SHAPE_CODE_SPHERE   0x53504845u // 'SPHE' in little-endian ASCII
#define TIN_SHAPE_CODE_POLYTOPE 0x50544F50u // 'PTOP' in little-endian ASCII

typedef uint32_t Tin_ShapeCode;
typedef struct Tin_Arbiter Tin_Arbiter;

typedef struct Tin_ShapeClass {
	Tin_ShapeCode code;
	Tin_Vec3 (*get_inv_inertia)(const void *shape);
	void (*get_aabb)(const void *shape, const Tin_Transform *transform, Tin_Vec3 *aabbMin, Tin_Vec3 *aabbMax);
	bool (*intersect)(const void *shape, const Tin_Transform *transform, const void *otherShape, const Tin_Transform *otherTransform, Tin_Arbiter *arbiter);
} Tin_ShapeClass;

typedef struct {
	const Tin_ShapeClass *vtable;
} Tin_Shape;

typedef struct {
	const Tin_ShapeClass *vtable;
	Tin_Scalar radius;
} Tin_Sphere;

void tin_shape_aabb(const Tin_Shape *shape, const Tin_Transform *transform, Tin_Vec3 *aabbMin, Tin_Vec3 *aabbMax);

extern const Tin_ShapeClass tin_shape_sphere;
extern const Tin_ShapeClass tin_shape_polytope;

/// @}

/// \defgroup body Rigid Bodies
/// @{

/// A numerical ID identifying a rigid body within a scene.
typedef uint32_t Tin_BodyID;

/// A rigid body.
typedef struct Tin_Body {
	/// The collision shape.
	const Tin_Shape *shape;
	/// The position, orientation, and scale in space.
	Tin_Transform transform;
	
	/// The linear velocity.
	Tin_Vec3      velocity;
	/**
	 * \brief The bodies' angular velocity.
	 *
	 * The direction of the vector represents the rotation axis,
	 * the magnitude represents the speed of the rotation in
	 * radians per second.
	 */
	Tin_Vec3      angularVelocity;
	
	/// One over the mass.
	Tin_Scalar    invMass;
	/// The inverse of the inertia tensor along the cardinal axes,
	/// assuming that the center of mass is at the origin of the rigid body.
	Tin_Scalar    invInertia[3*3];
	
	Tin_BodyID    islandID;
	int           restCounter;
	bool          islandStable;
	
	Tin_Vec3      aabbMin;
	Tin_Vec3      aabbMax;
} Tin_Body;

/// @}

/// \defgroup minkowski Minkowski Sum
/// @{

/// Implicit Minkowski Sum (difference) of two transformed polytopes.
typedef struct {
	const Tin_Polytope  *polytope1;
	const Tin_Transform *transform1;
	const Tin_Polytope  *polytope2;
	const Tin_Transform *transform2;
} Tin_Polysum;

Tin_Vec3 tin_polysum_support(const void *geometry, Tin_Vec3 dir);

/* === Minkowski Portal Refinement === :mpr: */

typedef struct Tin_Arbiter Tin_Arbiter;

/// A ray in 3D space.
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
bool tin_intersect(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Vec3 *normalOut);
void tin_fill_arbiter(
	const Tin_Transform *ta, const Tin_Transform *tb,
	int faceA, int faceB,
	Tin_Vec3 refNormal, Tin_Scalar refBase,
	Tin_Vec3 *manifold, int count,
	Tin_Arbiter *arbiter);
int tin_polytope_collide(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Arbiter *arbiter);

/// @}

/// \defgroup contact Contact Points
/// @{

#define TIN_MAX_CONTACTS 4

typedef struct {
	Tin_Scalar jacobian[12];
	Tin_Scalar invMassJacobian[12];
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

typedef struct Tin_Arbiter {
	int bodyID1;
	int bodyID2;
	int face1;
	int face2;
	Tin_Vec3 normal;
	Tin_Vec3 frictionDir;
	Tin_Vec3 orthoDir;
	int numContacts;
	int numPenetrating; // the penetrating contacts (separation < 0) come first
	Tin_Contact contacts[TIN_MAX_CONTACTS];
} Tin_Arbiter;

typedef struct Tin_Scene Tin_Scene;

void tin_arbiter_prestep(Tin_Scene *scene, Tin_Arbiter *arbiter, Tin_Scalar (*velocities)[6], Tin_Scalar invDt);

int tin_clip_manifold_against_plane(const Tin_Vec3 *points, int count, Tin_Vec3 normal, Tin_Scalar base, Tin_Vec3 *newPoints);
int tin_clip_manifolds(const Tin_Vec3 *pointsA, int countA, const Tin_Vec3 *pointsB, int countB, Tin_Vec3 perpendicular, Tin_Vec3 *newPoints);
int tin_incident_face(const Tin_Polytope *polytope, Tin_Vec3 dir);
int tin_reduce_manifold(Tin_Vec3 *points, int count);

/// @}

/// \defgroup pairtable Pair-Indexed Hashtable
/// @{

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
size_t tin_pairtable_index(const Tin_PairTable *table, size_t elemLow, size_t elemHigh);
bool tin_find_pair(const Tin_PairTable *table, size_t elemA, size_t elemB, void **payloadOut);
void tin_insert_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void *payload);
void tin_delete_pair(Tin_PairTable *table, size_t elemA, size_t elemB);

/// @}

/// \defgroup island Island Detection & Freezing
/// @{

typedef struct Tin_Scene Tin_Scene;

Tin_BodyID tin_island_find(Tin_Scene *scene, Tin_BodyID bodyID);
void tin_island_union(Tin_Scene *scene, Tin_BodyID bodyID1, Tin_BodyID bodyID2);
void tin_build_islands(Tin_Scene *scene);

/// @}

/// \defgroup scene Scenes
/// @{

typedef struct Tin_SweepPrune Tin_SweepPrune;

typedef struct Tin_Scene {
	Tin_Body  *bodyTable;
	bool      *bodyOccupied;
	Tin_BodyID freeBodyID;
	Tin_BodyID bodyTableCapac;

	Tin_Arbiter *arbiters;
	size_t numArbiters;
	size_t capArbiters;

	Tin_PairTable contactCache;
	Tin_Arbiter *oldArbiters;
	size_t numOldArbiters;
	size_t capOldArbiters;
	
	Tin_SweepPrune *sweepPrune;
} Tin_Scene;

void tin_scene_update(Tin_Scene *scene);
void tin_scene_prestep(Tin_Scene *scene, Tin_Scalar (*velocities)[6], Tin_Scalar invDt);
void tin_scene_step(Tin_Scene *scene, Tin_Scalar (*velocities)[6]);
void tin_integrate(Tin_Scene *scene, Tin_Scalar dt);
void tin_simulate(Tin_Scene *scene, Tin_Scalar dt, double (*gettime)(), double timings[6]);

/**
 * \brief Create a new rigid body.
 *
 * \param scene The scene to which the rigid body is added.
 * \param shape The shape of the rigid body.
 *              This pointer must remain valid for as long as the rigid body exists.
 *              It is not modified or freed by tincan.
 * \param invMass One over the mass of the rigid body. A value of zero means infinite mass.
 * \return A numerical ID identifying the body within the scene.
 */
Tin_BodyID tin_add_body(Tin_Scene *scene, const Tin_Shape *shape, Tin_Scalar invMass);

/**
 * \brief Delete a rigid body.
 *
 * \param scene  The scene from which the rigid body is to be deleted.
 * \param bodyID A numerical ID identifying the body within the scene.
 */
void       tin_delete_body(Tin_Scene *scene, Tin_BodyID bodyID);

/// @}

/// \defgroup broadphase Broadphase Collision Detection
/// @{

#define TIN_BLOOM_NUM_BITS (1024 * 8)

void tin_bloom_hash(uintptr_t key1, uintptr_t key2, unsigned hashes[3]);
void tin_bloom_insert(unsigned *bloom, uintptr_t key1, uintptr_t key2);
bool tin_bloom_lookup(const unsigned *bloom, uintptr_t key1, uintptr_t key2);

typedef struct {
	Tin_BodyID bodyID;
	Tin_Scalar value;
	bool isMin;
} Tin_SweepEvent;

typedef struct {
	size_t numEvents;
	size_t capEvents;
	Tin_SweepEvent *events;
} Tin_SweepAxis;

/// Incremental Sweep-and-Prune broadphase collision detection.
typedef struct Tin_SweepPrune {
	Tin_Scene *scene;
	Tin_SweepAxis axes[3];
} Tin_SweepPrune;

void tin_create_sweep_prune(Tin_SweepPrune *sap, Tin_Scene *scene);
void tin_destroy_sweep_prune(Tin_SweepPrune *sap);
void tin_sweep_prune_add_body(Tin_SweepPrune *sap, Tin_BodyID bodyID);
void tin_sweep_prune_delete_body(Tin_SweepPrune *sap, Tin_BodyID bodyID);
void tin_sweep_prune_update(Tin_SweepPrune *sap);
void tin_sweep_prune_axis(Tin_SweepPrune *sap, int axisIdx, const unsigned *bloomFilter,
	Tin_BodyID **collidersOut, size_t *numCollidersOut);
void tin_sweep_prune(Tin_SweepPrune *sap,
	Tin_BodyID **collidersOut, size_t *numCollidersOut);

/// @}

#endif
