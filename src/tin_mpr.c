#include "tincan.h"

#include <stdio.h>
#include <string.h>
#include <tgmath.h>

int
tin_construct_portal(const void *geometry, Tin_SupportFunc support,
	const Tin_Ray *r, Tin_Portal *p)
{
	/* find the very first point */
	p->a = support(geometry, r->dir);

	/* find the second point */
	{
		Tin_Vec3 oa = tin_sub_v3(p->a, r->origin);
		Tin_Vec3 dir = tin_gram_schmidt(oa, r->dir);
		if (tin_dot_v3(dir, dir) == 0.0f) {
			dir = tin_gram_schmidt(oa, (Tin_Vec3) {{ 1.0f, 0.0f, 0.0f }});
			if (tin_dot_v3(dir, dir) == 0.0f) {
				dir = (Tin_Vec3) {{ 0.0f, 0.0f, 1.0f }};
			}
		}
		dir = tin_normalize_v3(dir);
		p->b = support(geometry, dir);
	}

	for (int it = 0;; it++) {
		if (it >= 100) {
			fprintf(stderr, "Portal construction took too many iterations.\n");
			return 0;
		}

		/*  */
		Tin_Vec3 oa = tin_sub_v3(p->a, r->origin);
		Tin_Vec3 ob = tin_sub_v3(p->b, r->origin);
		Tin_Vec3 dir = tin_cross_v3(oa, ob);
		if (tin_dot_v3(dir, r->dir) < 0.0f) {
			dir = tin_neg_v3(dir);
		}
		p->c = support(geometry, dir);
		Tin_Scalar score = tin_dot_v3(dir, p->c);
		if (score <= tin_dot_v3(dir, p->a) || score <= tin_dot_v3(dir, p->b)) {
			return 0;
		}
		Tin_Vec3 oc = tin_sub_v3(p->c, r->origin);
		
		/*  */
		Tin_Vec3 an = tin_cross_v3(ob, oc);
		if (tin_dot_v3(an, oa) * tin_dot_v3(an, r->dir) < 0.0f) {
			p->a = p->c;
			continue;
		}
		
		/*  */
		Tin_Vec3 bn = tin_cross_v3(oc, oa);
		if (tin_dot_v3(bn, ob) * tin_dot_v3(bn, r->dir) < 0.0f) {
			p->b = p->c;
			continue;
		}
		
		return 1;
	}
}

void
tin_refine_portal(const void *geometry, Tin_SupportFunc support,
	const Tin_Ray *r, Tin_Portal *p)
{
	Tin_Vec3 ab = tin_sub_v3(p->b, p->a);
	Tin_Vec3 ac = tin_sub_v3(p->c, p->a);
	p->normal = tin_cross_v3(ab, ac);
	if (tin_dot_v3(p->normal, r->dir) < 0.0f) {
		Tin_Vec3 temp;
		temp = p->a;
		p->a = p->b;
		p->b = temp;
	}

	for (int it = 0;; it++) {
		Tin_Vec3 ab = tin_sub_v3(p->b, p->a);
		Tin_Vec3 ac = tin_sub_v3(p->c, p->a);
		p->normal = tin_cross_v3(ab, ac);

		if (tin_dot_v3(p->normal, p->normal) == 0.0f) {
			fprintf(stderr, "portal collapsed. (it=%d)\n", it);
			break;
		}

		Tin_Vec3 s = support(geometry, p->normal);
		
		{
			Tin_Scalar sScore = tin_dot_v3(s, p->normal);
			if (sScore <= 0.0f) break;
			Tin_Scalar aScore = tin_dot_v3(p->a, p->normal);
			if (sScore <= aScore * 1.01f) break;
			if (sScore - aScore <= 1e-6f) break;
		}
		
		if (it >= 100) {
			fprintf(stderr, "MPR refine_portal() took too many iterations.\n");
			fprintf(stderr, "\tnormal length squared: %f\n", tin_dot_v3(p->normal, p->normal));
			Tin_Scalar aScore = tin_dot_v3(p->a, p->normal);
			Tin_Scalar sScore = tin_dot_v3(s, p->normal);
			fprintf(stderr, "\taScore: %e\n", aScore);
			fprintf(stderr, "\tsScore: %e\n", sScore);
			break;
		}

		Tin_Vec3 os = tin_sub_v3(s, r->origin);
		Tin_Vec3 d_x_os = tin_cross_v3(r->dir, os);
		
		/* if <D , (OS x OA)> > 0 */
		Tin_Vec3 oa = tin_sub_v3(p->a, r->origin);
		if (tin_dot_v3(oa, d_x_os) > 0.0f) {
			/* if <D , (OS x OB)> > 0 */
			Tin_Vec3 ob = tin_sub_v3(p->b, r->origin);
			if (tin_dot_v3(ob, d_x_os) > 0.0f) {
				p->a = s; /* SBC */
			} else {
				p->c = s; /* ABS */
			}
		} else {
			/* if <D , (OS x OC)> > 0 */
			Tin_Vec3 oc = tin_sub_v3(p->c, r->origin);
			if (tin_dot_v3(oc, d_x_os) > 0.0f) {
				p->b = s; /* ASC */
			} else {
				p->a = s; /* SBC */
			}
		}
	}
}

bool
tin_intersect(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Vec3 *normalOut)
{
	Tin_Polysum ps = { pa, ta, pb, tb };

	Tin_Ray r;
	r.origin = tin_sub_v3(ta->translation, tb->translation);
	r.dir    = tin_neg_v3(r.origin);
	Tin_Scalar norm = sqrt(tin_dot_v3(r.dir, r.dir));
	if (norm == 0.0) {
		/* FIXME */
		return false;
	}
	r.dir = tin_normalize_v3(r.dir);

	Tin_Portal p;
	if (!tin_construct_portal(&ps, tin_polysum_support, &r, &p)) {
		return false;
	}
	tin_refine_portal(&ps, tin_polysum_support, &r, &p);
	p.normal = tin_normalize_v3(p.normal);
	if (tin_dot_v3(p.normal, p.a) <= 0.0) {
		return false;
	}

	if (normalOut == NULL) {
		return true;
	}

	for (int it = 0; it < 4; it++) {
		Tin_Ray nr;
		nr.dir    = p.normal;
		nr.origin = TIN_VEC3(0.0, 0.0, 0.0);
		Tin_Portal np;
		if (!tin_construct_portal(&ps, tin_polysum_support, &nr, &np)) {
			break;
		}
		tin_refine_portal(&ps, tin_polysum_support, &nr, &np);
		np.normal = tin_normalize_v3(np.normal);

		Tin_Scalar proj = tin_dot_v3(p.normal, np.normal);
		p = np; /* lol */
		r = nr;
		if (proj >= 0.99) break;
	}

	*normalOut = p.normal;
	return true;
}

void
tin_fill_arbiter(
	const Tin_Transform *ta, const Tin_Transform *tb,
	int faceA, int faceB,
	Tin_Vec3 refNormal, Tin_Scalar refBase,
	Tin_Vec3 *manifold, int count,
	Tin_Arbiter *arbiter)
{
	// Fill out the Arbiter structure with contact point information.
	arbiter->normal = refNormal;
	arbiter->face1 = faceA;
	arbiter->face2 = faceB;
	memset(arbiter->contacts, 0, count * sizeof *arbiter->contacts);
	for (int i = 0; i < count; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];
		Tin_Scalar projDist = tin_dot_v3(refNormal, manifold[i]) - refBase;
		Tin_Vec3 projPoint = tin_saxpy_v3(-projDist, refNormal, manifold[i]);
		contact->posFrom1 = tin_sub_v3(projPoint, ta->translation);
		contact->posFrom2 = tin_sub_v3(manifold[i], tb->translation);
		Tin_Vec3 p1 = projPoint;
		Tin_Vec3 p2 = manifold[i];
		contact->separation = tin_dot_v3(refNormal, tin_sub_v3(p2, p1));
	}

	// Swap around contact points so that the ones that penetrate are at the front of the array.
	for (int n = count; n > 1; n--) {
		for (int i = 0; i < n - 1; i++) {
			if (arbiter->contacts[i].separation >= 0.0 && !(arbiter->contacts[i+1].separation >= 0.0)) {
				Tin_Contact temp = arbiter->contacts[i];
				arbiter->contacts[i] = arbiter->contacts[i+1];
				arbiter->contacts[i+1] = temp;
			}
		}
	}

	// Find out how many penetrating contact points there are.
	int numPenetrating = 0;
	for (int i = 0; i < count; i++) {
		if (!(arbiter->contacts[i].separation >= 0.0)) {
			numPenetrating++;
		}
	}
	arbiter->numPenetrating = numPenetrating;
}

