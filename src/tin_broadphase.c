#include "tincan.h"

#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#define SET_BIT(u,i) ((u)[(i) / (8*sizeof*(u))] |= 1U << ((i) % (8*sizeof*(u))))
#define GET_BIT(u,i) (((u)[(i) / (8*sizeof*(u))] >> ((i) % (8*sizeof*(u)))) & 1U)

void
tin_bloom_hash(uintptr_t key1, uintptr_t key2, unsigned hashes[3])
{
	uint64_t x;
	if (key1 < key2) {
		x = (key1 << 32) ^ key2;
	} else {
		x = (key2 << 32) ^ key1;
	}

	/* Some 64-bit bijective hash function, taken from:
	 * https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
	 */
	x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
	x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
	x = x ^ (x >> 31);

	hashes[0] = x & 0x1FFFFF;
	hashes[1] = (x >> 21) & 0x1FFFFF;
	hashes[2] = (x >> 42) & 0x1FFFFF;
}

void
tin_bloom_insert(unsigned *bloom, uintptr_t key1, uintptr_t key2)
{
	unsigned idx[3];
	tin_bloom_hash(key1, key2, idx);
	idx[0] %= TIN_BLOOM_NUM_BITS;
	idx[1] %= TIN_BLOOM_NUM_BITS;
	idx[2] %= TIN_BLOOM_NUM_BITS;
	SET_BIT(bloom, idx[0]);
	SET_BIT(bloom, idx[1]);
	SET_BIT(bloom, idx[2]);
}

bool
tin_bloom_lookup(const unsigned *bloom, uintptr_t key1, uintptr_t key2)
{
	unsigned idx[3];
	tin_bloom_hash(key1, key2, idx);
	idx[0] %= TIN_BLOOM_NUM_BITS;
	idx[1] %= TIN_BLOOM_NUM_BITS;
	idx[2] %= TIN_BLOOM_NUM_BITS;
	return GET_BIT(bloom, idx[0]) & GET_BIT(bloom, idx[1]) & GET_BIT(bloom, idx[2]);
}

void
tin_create_sweep_prune(Tin_SweepPrune *sap, Tin_Scene *scene)
{
	sap->scene = scene;
	for (int a = 0; a < 3; a++) {
		sap->axes[a].numEvents = 0;
		sap->axes[a].events = NULL;
	}
}

void
tin_destroy_sweep_prune(Tin_SweepPrune *sap)
{
	for (int a = 0; a < 3; a++) {
		free(sap->axes[a].events);
	}
}

void
tin_sweep_prune_add_body(Tin_SweepPrune *sap, Tin_BodyID bodyID)
{
	for (int a = 0; a < 3; a++) {
		Tin_SweepAxis *axis = &sap->axes[a];
		if (axis->numEvents + 2 > axis->capEvents) {
			axis->capEvents *= 2;
			if (!axis->capEvents) axis->capEvents = 16;
			axis->events = realloc(axis->events, axis->capEvents * sizeof *axis->events);
		}
		axis->events[axis->numEvents++] = (Tin_SweepEvent){
			.bodyID =  bodyID,
			.value  = -INFINITY,
			.isMin  =  true
		};
		axis->events[axis->numEvents++] = (Tin_SweepEvent){
			.bodyID =  bodyID,
			.value  =  INFINITY,
			.isMin  =  false
		};
	}
}

void
tin_sweep_prune_delete_body(Tin_SweepPrune *sap, Tin_BodyID bodyID)
{
	for (int a = 0; a < 3; a++) {
		Tin_SweepAxis *axis = &sap->axes[a];
		for (size_t e = 0; e < axis->numEvents; e++) {
			if (axis->events[e].bodyID == bodyID) {
				axis->events[e] = axis->events[--axis->numEvents];
				// If event count is even again, we must have deleted
				// both events, so we can exit early.
				if ((axis->numEvents & 1) == 0) break;
			}
		}
	}
}

void
tin_sweep_prune_update(Tin_SweepPrune *sap)
{
	/* Update event positions (body extents) */
	for (int a = 0; a < 3; a++) {
		Tin_SweepAxis *axis = &sap->axes[a];
		for (size_t e = 0; e < axis->numEvents; e++) {
			Tin_SweepEvent *event = &axis->events[e];
			const Tin_Body *body = &sap->scene->bodyTable[event->bodyID];
			Tin_Vec3 origin = body->transform.translation;
			Tin_Scalar radius = body->shape->radius * body->transform.scale;
			event->value = origin.c[a] + (event->isMin ? -radius : radius);
		}
	}

	/* Adapative, in-place event sorting (insertion sort) */
	for (int a = 0; a < 3; a++) {
		Tin_SweepAxis *axis = &sap->axes[a];
		for (size_t j, i = 1; i < axis->numEvents; i++) {
			Tin_SweepEvent temp = axis->events[i];
			for (j = i; j > 0 && axis->events[j-1].value > temp.value; j--) {
				axis->events[j] = axis->events[j-1];
			}
			axis->events[j] = temp;
		}
	}
}

void
tin_sweep_prune_axis(Tin_SweepPrune *sap, int axisIdx, const unsigned *bloomFilter,
	Tin_BodyID **collidersOut, size_t *numCollidersOut)
{
	const Tin_SweepAxis *axis = &sap->axes[axisIdx];

	Tin_BodyID *colliders = NULL;
	size_t numColliders = 0;
	size_t capColliders = 0;

	Tin_BodyID *active = calloc(axis->numEvents / 2, sizeof *active);
	size_t numActive = 0;

	for (size_t e = 0; e < axis->numEvents; e++) {
		const Tin_SweepEvent *event = &axis->events[e];
		if (event->isMin) {
			for (size_t a = 0; a < numActive; a++) {
				if (numColliders + 2 > capColliders) {
					capColliders *= 2;
					if (!capColliders) capColliders = 16;
					colliders = realloc(colliders, capColliders * sizeof *colliders);
				}
				if (!bloomFilter || tin_bloom_lookup(bloomFilter, active[a], event->bodyID)) {
					colliders[numColliders++] = active[a];
					colliders[numColliders++] = event->bodyID;
				}
			}
			active[numActive++] = event->bodyID;
		} else {
			size_t a;
			for (a = 0; a < numActive; a++) {
				if (active[a] == event->bodyID) break;
			}
			active[a] = active[--numActive];
		}
	}

	free(active);

	*collidersOut = colliders;
	*numCollidersOut = numColliders;
}

void
tin_sweep_prune(Tin_SweepPrune *sap,
	Tin_BodyID **collidersOut, size_t *numCollidersOut)
{
	Tin_BodyID *colliders;
	size_t numColliders;

	tin_sweep_prune_axis(sap, 0, NULL, &colliders, &numColliders);

	unsigned *bloomFilter = calloc(1, TIN_BLOOM_NUM_BITS / 8);
	for (size_t i = 0; i < numColliders; i += 2) {
		tin_bloom_insert(bloomFilter, colliders[i], colliders[i+1]);
	}
	free(colliders);

	tin_sweep_prune_axis(sap, 1, bloomFilter, &colliders, &numColliders);

	memset(bloomFilter, 0, TIN_BLOOM_NUM_BITS / 8);
	for (size_t i = 0; i < numColliders; i += 2) {
		tin_bloom_insert(bloomFilter, colliders[i], colliders[i+1]);
	}
	free(colliders);

	tin_sweep_prune_axis(sap, 2, bloomFilter, &colliders, &numColliders);

	free(bloomFilter);

	/* Eliminate false positives */
	size_t w = 0;
	for (size_t i = 0; i < numColliders; i += 2) {
		const Tin_Body *body1 = &sap->scene->bodyTable[colliders[i+0]];
		const Tin_Body *body2 = &sap->scene->bodyTable[colliders[i+1]];
		Tin_Vec3 diff = tin_sub_v3(body2->transform.translation, body1->transform.translation);
		Tin_Scalar radii =
			body2->shape->radius * body2->transform.scale +
			body1->shape->radius * body1->transform.scale;
		if (tin_dot_v3(diff, diff) <= radii * radii) {
			colliders[w++] = colliders[i+0];
			colliders[w++] = colliders[i+1];
		}
	}
	numColliders = w;

	*collidersOut = colliders;
	*numCollidersOut = numColliders;
}
