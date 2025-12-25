#include "tincan.h"

#include <stdlib.h>
#include <string.h>

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
	union {
		uint64_t q[2];
		unsigned char b[16];
	} u;
	u.q[0] = elemLow;
	u.q[1] = elemHigh;

	/* FNV-1a hashing */
	uint32_t hash = 2166136261u;
	for (int i = 0; i < 16; i++) {
		hash ^= u.b[i];
		hash *= 16777619u;
	}
	return hash;
}

size_t
tin_fold_hash(size_t capac, uint32_t hash)
{
	// TODO xor folding or similar
	uint32_t mask = capac - 1;
	return hash & mask;
}

size_t
tin_pairtable_index(const Tin_PairTable *table, size_t elemLow, size_t elemHigh)
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
tin_find_pair(const Tin_PairTable *table, size_t elemA, size_t elemB, void **payloadOut)
{
	tin_order_pair(&elemA, &elemB);
	size_t index = tin_pairtable_index(table, elemA, elemB);
	if (!table->slots[index].occupied) {
		return false;
	}
	if (payloadOut) {
		*payloadOut = table->slots[index].payload;
	}
	return true;
}

void
tin_insert_pair(Tin_PairTable *table, size_t elemA, size_t elemB, void *payload)
{
	if (2 * table->count >= table->capac) {
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

