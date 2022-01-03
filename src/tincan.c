#include "tincan.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

tin_vec3
tin_neg_v3(tin_vec3 x)
{
	return (tin_vec3) {{ -x.c[0], -x.c[1], -x.c[2] }};
}

tin_vec3
tin_add_v3(tin_vec3 a, tin_vec3 b)
{
	return (tin_vec3) {{ a.c[0]+b.c[0], a.c[1]+b.c[1], a.c[2]+b.c[2] }};
}

tin_vec3
tin_sub_v3(tin_vec3 a, tin_vec3 b)
{
	return (tin_vec3) {{ a.c[0]-b.c[0], a.c[1]-b.c[1], a.c[2]-b.c[2] }};
}

tin_vec3
tin_scale_v3(tin_scalar a, tin_vec3 x)
{
	return (tin_vec3) {{ a*x.c[0], a*x.c[1], a*x.c[2] }};
}

tin_vec3
tin_saxpy_v3(tin_scalar a, tin_vec3 x, tin_vec3 y)
{
	return tin_add_v3(tin_scale_v3(a, x), y);
}

tin_vec3
tin_cross_v3(tin_vec3 a, tin_vec3 b)
{
	tin_vec3 c;
	c.c[0] = a.c[1] * b.c[2] - a.c[2] * b.c[1];
	c.c[1] = a.c[2] * b.c[0] - a.c[0] * b.c[2];
	c.c[2] = a.c[0] * b.c[1] - a.c[1] * b.c[0];
	return c;
}

tin_scalar
tin_dot_v3(tin_vec3 a, tin_vec3 b)
{
	return a.c[0]*b.c[0] + a.c[1]*b.c[1] + a.c[2]*b.c[2];
}

tin_vec3
tin_normalize_v3(tin_vec3 v)
{
	tin_scalar norm = sqrtf(tin_dot_v3(v, v));
	if (norm > 0.0f) {
		return tin_scale_v3(1.0f / norm, v);
	} else {
		return v;
	}
}

tin_vec3
tin_hadamard_v3(tin_vec3 a, tin_vec3 b)
{
	return (tin_vec3) {{ a.c[0]*b.c[0], a.c[1]*b.c[1], a.c[2]*b.c[2] }};
}

tin_quat
tin_make_qt(tin_vec3 axis, tin_scalar angle)
{
	tin_scalar h = angle / 2.0f;
	return (tin_quat) { tin_scale_v3(sinf(h), axis), cosf(h) };
}

tin_quat
tin_mul_qt(tin_quat a, tin_quat b)
{
	tin_quat d;
	d.v = tin_cross_v3(a.v, b.v);
	d.v = tin_saxpy_v3(a.s, b.v, d.v);
	d.v = tin_saxpy_v3(b.s, a.v, d.v);
	d.s = a.s * b.s - tin_dot_v3(a.v, b.v);
	return d;
}

tin_vec3
tin_apply_qt(tin_quat q, tin_vec3 v)
{
	tin_vec3 t = tin_cross_v3(q.v, v);
	tin_vec3 d = tin_cross_v3(q.v, t);
	d = tin_saxpy_v3(q.s, t, d);
	d = tin_saxpy_v3(2.0f, d, v);
	return d;
}

tin_quat
tin_conjugate_qt(tin_quat q)
{
	return (tin_quat) { tin_neg_v3(q.v), q.s };
}

tin_vec3
tin_gram_schmidt(tin_vec3 fixed, tin_vec3 var)
{
	tin_scalar fscale = tin_dot_v3(fixed, fixed);
	if (fscale) {
		return tin_saxpy_v3(-tin_dot_v3(fixed, var) / fscale, fixed, var);
	} else {
		return var;
	}
}

tin_scalar
tin_prlgram_area(tin_vec3 e1, tin_vec3 e2)
{
	tin_vec3 perp = tin_gram_schmidt(e1, e2);
	return sqrtf(tin_dot_v3(e1, e1) * tin_dot_v3(perp, perp));
}

tin_vec3
tin_fwtrf_point(const tin_transform *transform, tin_vec3 vec)
{
	vec = tin_apply_qt(transform->rotation,    vec);
	vec = tin_scale_v3(transform->scale,       vec);
	vec = tin_add_v3  (transform->translation, vec);
	return vec;
}

tin_vec3
tin_fwtrf_dir(const tin_transform *transform, tin_vec3 vec)
{
	return tin_apply_qt(transform->rotation, vec);
}

tin_vec3
tin_bwtrf_point(const tin_transform *transform, tin_vec3 vec)
{
	vec = tin_sub_v3  (vec, transform->translation);
	vec = tin_scale_v3(1.0f / transform->scale, vec);
	vec = tin_apply_qt(tin_conjugate_qt(transform->rotation), vec);
	return vec;
}

tin_vec3
tin_bwtrf_dir(const tin_transform *transform, tin_vec3 vec)
{
	return tin_apply_qt(tin_conjugate_qt(transform->rotation), vec);
}

tin_vec3
tin_polytope_support(const tin_polytope *p, tin_vec3 dir)
{
	tin_scalar score, best_score = -INFINITY;
	int i, best_i = -1;
	for (i = 0; i < p->num_vertices; i++) {
		score = tin_dot_v3(p->vertices[i], dir);
		if (score > best_score) {
			best_score = score;
			best_i     = i;
		}
	}
	return p->vertices[best_i];
}

void
tin_polysum_support(const tin_polysum *s, tin_vec3 dir, tin_pspoint *sup)
{
	tin_vec3 former_dir = tin_bwtrf_dir(s->former_transform, dir);
	sup->rel_former = tin_polytope_support(s->former_polytope, former_dir);
	tin_vec3 former_abs = tin_fwtrf_point(s->former_transform, sup->rel_former);

	tin_vec3 latter_dir = tin_bwtrf_dir(s->latter_transform, tin_neg_v3(dir));
	sup->rel_latter = tin_polytope_support(s->latter_polytope, latter_dir);
	tin_vec3 latter_abs = tin_fwtrf_point(s->latter_transform, sup->rel_latter);

	sup->abs = tin_sub_v3(former_abs, latter_abs);
}

int
tin_construct_portal(const tin_polysum *ps, const tin_ray *r, tin_portal *p)
{
	/* find the very first point */
	tin_polysum_support(ps, r->dir, &p->a);

	/* find the second point */
	{
		tin_vec3 oa = tin_sub_v3(p->a.abs, r->origin);
		tin_vec3 dir = tin_gram_schmidt(oa, r->dir);
		if (tin_dot_v3(dir, dir) == 0.0f) {
			dir = tin_gram_schmidt(oa, (tin_vec3) {{ 1.0f, 0.0f, 0.0f }});
			if (tin_dot_v3(dir, dir) == 0.0f) {
				dir = (tin_vec3) {{ 0.0f, 0.0f, 1.0f }};
			}
		}
		dir = tin_normalize_v3(dir);
		tin_polysum_support(ps, dir, &p->b);
	}

	for (int it = 0;; it++) {
		if (it >= 100) {
			fprintf(stderr, "Portal construction took too many iterations.\n");
			return 0;
		}

		/*  */
		tin_vec3 oa = tin_sub_v3(p->a.abs, r->origin);
		tin_vec3 ob = tin_sub_v3(p->b.abs, r->origin);
		tin_vec3 dir = tin_cross_v3(oa, ob);
		if (tin_dot_v3(dir, r->dir) < 0.0f) {
			dir = tin_neg_v3(dir);
		}
		tin_polysum_support(ps, dir, &p->c);
		tin_scalar score = tin_dot_v3(dir, p->c.abs);
		if (score <= tin_dot_v3(dir, p->a.abs) || score <= tin_dot_v3(dir, p->b.abs)) {
			return 0;
		}
		tin_vec3 oc = tin_sub_v3(p->c.abs, r->origin);
		
		/*  */
		tin_vec3 an = tin_cross_v3(ob, oc);
		if (tin_dot_v3(an, oa) * tin_dot_v3(an, r->dir) < 0.0f) {
			p->a = p->c;
			continue;
		}
		
		/*  */
		tin_vec3 bn = tin_cross_v3(oc, oa);
		if (tin_dot_v3(bn, ob) * tin_dot_v3(bn, r->dir) < 0.0f) {
			p->b = p->c;
			continue;
		}
		
		return 1;
	}
}

void
tin_refine_portal(const tin_polysum *ps, const tin_ray *r, tin_portal *p)
{
	tin_vec3 ab = tin_sub_v3(p->b.abs, p->a.abs);
	tin_vec3 ac = tin_sub_v3(p->c.abs, p->a.abs);
	p->normal = tin_cross_v3(ab, ac);
	if (tin_dot_v3(p->normal, r->dir) < 0.0f) {
		tin_pspoint temp;
		temp = p->a;
		p->a = p->b;
		p->b = temp;
	}

	for (int it = 0;; it++) {
		tin_vec3 ab = tin_sub_v3(p->b.abs, p->a.abs);
		tin_vec3 ac = tin_sub_v3(p->c.abs, p->a.abs);
		p->normal = tin_cross_v3(ab, ac);

		if (tin_dot_v3(p->normal, p->normal) == 0.0f) {
			fprintf(stderr, "portal collapsed. (it=%d)\n", it);
			break;
		}

		if (it >= 100) {
			fprintf(stderr, "MPR refine_portal() took too many iterations.\n");
			break;
		}

		tin_pspoint s;
		tin_polysum_support(ps, p->normal, &s);
		
		{
			tin_vec3 as = tin_sub_v3(s.abs, p->a.abs);
			tin_vec3 bs = tin_sub_v3(s.abs, p->b.abs);
			tin_vec3 cs = tin_sub_v3(s.abs, p->c.abs);
			if (tin_dot_v3(as, p->normal) <= 0.0f) {
				break;
			}
			if (tin_dot_v3(bs, p->normal) <= 0.0f) {
				break;
			}
			if (tin_dot_v3(cs, p->normal) <= 0.0f) {
				break;
			}
		}
		
		tin_vec3 os = tin_sub_v3(s.abs, r->origin);
		tin_vec3 d_x_os = tin_cross_v3(r->dir, os);
		
		/* if <D , (OS x OA)> > 0 */
		tin_vec3 oa = tin_sub_v3(p->a.abs, r->origin);
		if (tin_dot_v3(oa, d_x_os) > 0.0f) {
			/* if <D , (OS x OB)> > 0 */
			tin_vec3 ob = tin_sub_v3(p->b.abs, r->origin);
			if (tin_dot_v3(ob, d_x_os) > 0.0f) {
				p->a = s; /* SBC */
			} else {
				p->c = s; /* ABS */
			}
		} else {
			/* if <D , (OS x OC)> > 0 */
			tin_vec3 oc = tin_sub_v3(p->c.abs, r->origin);
			if (tin_dot_v3(oc, d_x_os) > 0.0f) {
				p->b = s; /* ASC */
			} else {
				p->a = s; /* SBC */
			}
		}
	}
}

void
tin_calculate_contact(const tin_ray *r, const tin_portal *p, tin_contact *contact)
{
	contact->normal = tin_normalize_v3(p->normal);

	tin_vec3 ab = tin_sub_v3(p->b.abs, p->a.abs);
	tin_vec3 ac = tin_sub_v3(p->c.abs, p->a.abs);
	tin_scalar area_total = tin_prlgram_area(ab, ac);
	if (area_total == 0.0f) {
		contact->rel1 = p->a.rel_former;
		contact->rel2 = p->a.rel_latter;
		fprintf(stderr, "Calculating contact point on collapsed portal.\n");
		return;
	}
	tin_scalar proj = tin_dot_v3(p->normal, r->dir);
	if (proj == 0.0f) {
		contact->rel1 = p->a.rel_former;
		contact->rel2 = p->a.rel_latter;
		fprintf(stderr, "Portal and contact normal have become parallel.\n");
		return;
	}
	tin_scalar t = tin_dot_v3(p->normal, tin_sub_v3(p->a.abs, r->origin)) / proj;
	tin_vec3 q = tin_saxpy_v3(t, r->dir, r->origin);
	tin_vec3 aq = tin_sub_v3(q, p->a.abs);
	tin_scalar area_b = tin_prlgram_area(aq, ac);
	tin_scalar area_c = tin_prlgram_area(aq, ab);
	tin_scalar beta = area_b / area_total;
	tin_scalar gamma = area_c / area_total;
	tin_scalar alpha = 1.0f - beta - gamma;
	
	contact->rel1 = tin_scale_v3(alpha, p->a.rel_former);
	contact->rel1 = tin_saxpy_v3(beta,  p->b.rel_former, contact->rel1);
	contact->rel1 = tin_saxpy_v3(gamma, p->c.rel_former, contact->rel1);

	contact->rel2 = tin_scale_v3(alpha, p->a.rel_latter);
	contact->rel2 = tin_saxpy_v3(beta,  p->b.rel_latter, contact->rel2);
	contact->rel2 = tin_saxpy_v3(gamma, p->c.rel_latter, contact->rel2);
}

int
tin_polytope_collide(
	const tin_polytope *pa, const tin_transform *ta,
	const tin_polytope *pb, const tin_transform *tb,
	tin_contact *contact)
{
	tin_polysum ps = { pa, ta, pb, tb };

	tin_ray r;
	r.origin = tin_sub_v3(ta->translation, tb->translation);
	r.dir    = tin_neg_v3(r.origin);
	tin_scalar norm = sqrtf(tin_dot_v3(r.dir, r.dir));
	if (norm == 0.0f) {
		/* FIXME */
		return 1;
	}
	r.dir = tin_normalize_v3(r.dir);

	tin_portal p;
	if (!tin_construct_portal(&ps, &r, &p)) {
		return 0;
	}
	tin_refine_portal(&ps, &r, &p);
	p.normal = tin_normalize_v3(p.normal);
	if (tin_dot_v3(p.normal, p.a.abs) <= 0.0f) {
		return 0;
	}
	p.normal = tin_normalize_v3(p.normal);

	for (int it = 0; it < 4; it++) {
		tin_ray nr;
		nr.dir    = p.normal;
		nr.origin = (tin_vec3) {{ 0.0f, 0.0f, 0.0f }};
		tin_portal np;
		if (!tin_construct_portal(&ps, &nr, &np)) {
			break;
		}
		tin_refine_portal(&ps, &nr, &np);
		np.normal = tin_normalize_v3(np.normal);

		tin_scalar proj = tin_dot_v3(p.normal, np.normal);
		p = np; /* unintentional, i swear */
		r = nr;
		if (proj >= 0.99f) break;
	}

	tin_calculate_contact(&r, &p, contact);
	return 1;
}

static tin_vec3
tin_solve_inertia(const tin_body *b, tin_vec3 vec)
{
	vec = tin_bwtrf_dir(&b->transform, vec);
	vec = tin_hadamard_v3(b->inv_inertia, vec);
	vec = tin_fwtrf_dir(&b->transform, vec);
	return vec;
}

static tin_scalar
tin_length_v3(tin_vec3 v)
{
	return sqrtf(tin_dot_v3(v, v));
}

void
tin_arbiter_update(tin_arbiter *a)
{
	const tin_scalar max_separation = 0.1f;
	const tin_scalar max_stretch = 0.3f;

	tin_body *b1 = a->body1;
	tin_body *b2 = a->body2;

	for (int i = 0; i < a->num_contacts; i++) {
		tin_contact *c = &a->contacts[i];

		tin_vec3 p1 = tin_fwtrf_point(&b1->transform, c->rel1);
		tin_vec3 p2 = tin_fwtrf_point(&b2->transform, c->rel2);

		tin_scalar separation = tin_dot_v3(c->normal, tin_sub_v3(p2, p1));

		if (separation > max_separation) {
			*c = a->contacts[--a->num_contacts];
			i--;
			continue;
		}

		tin_scalar stretch = tin_length_v3(tin_gram_schmidt(c->normal, tin_sub_v3(p1, p2)));
		if (stretch > max_stretch) {
			*c = a->contacts[--a->num_contacts];
			i--;
			continue;
		}
	}

}

void
tin_arbiter_add_contact(tin_arbiter *a, tin_contact contact)
{
	tin_body *b1 = a->body1;
	tin_body *b2 = a->body2;

	tin_vec3 p1, p2;
	p1 = tin_fwtrf_point(&b1->transform, contact.rel1);
	p2 = tin_fwtrf_point(&b2->transform, contact.rel2);
	contact.base_stretch = tin_length_v3(tin_gram_schmidt(contact.normal, tin_sub_v3(p1, p2)));

	if (a->num_contacts < TIN_MAX_CONTACTS) {
		a->contacts[a->num_contacts++] = contact;
		return;
	}

	tin_vec3 new_pos = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
	int best_i = -1;
	tin_scalar best_dist = INFINITY;
	for (int i = 0; i < a->num_contacts; i++) {
		tin_contact *c = &a->contacts[i];
		p1 = tin_fwtrf_point(&b1->transform, c->rel1);
		p2 = tin_fwtrf_point(&b2->transform, c->rel2);
		tin_vec3 old_pos = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
		tin_vec3 diff = tin_sub_v3(old_pos, new_pos);
		tin_scalar dist = tin_dot_v3(diff, diff);
		if (dist < best_dist) {
			best_i = i;
			best_dist = dist;
		}
	}
	a->contacts[best_i] = contact;
}

void
tin_arbiter_prestep(tin_arbiter *a, tin_scalar inv_dt)
{
	const tin_scalar allowed_penetration = 0.01f;
	const tin_scalar bias_factor = 0.2f;

	tin_body *b1 = a->body1;
	tin_body *b2 = a->body2;

	for (int i = 0; i < a->num_contacts; i++) {
		tin_contact *c = &a->contacts[i];

		tin_vec3 p1 = tin_fwtrf_point(&b1->transform, c->rel1);
		tin_vec3 p2 = tin_fwtrf_point(&b2->transform, c->rel2);

		c->position   = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
		c->separation = tin_dot_v3(c->normal, tin_sub_v3(p2, p1));

		tin_vec3 r1 = tin_sub_v3(c->position, b1->transform.translation);
		tin_vec3 r2 = tin_sub_v3(c->position, b2->transform.translation);

		/* Precompute normal mass and bias */
		c->normal_mass = b1->inv_mass + b2->inv_mass + tin_dot_v3(c->normal, tin_add_v3(
			tin_cross_v3(tin_solve_inertia(b1, tin_cross_v3(r1, c->normal)), r1),
			tin_cross_v3(tin_solve_inertia(b2, tin_cross_v3(r2, c->normal)), r2)));

		c->bias = -bias_factor * inv_dt * MIN(0.0f, c->separation + allowed_penetration);
	}
}

void
tin_arbiter_apply_impulse(tin_arbiter *a)
{
	tin_body *b1 = a->body1;
	tin_body *b2 = a->body2;

	for (int i = 0; i < a->num_contacts; i++) {
		tin_contact *c = &a->contacts[i];

		tin_vec3 r1 = tin_sub_v3(c->position, b1->transform.translation);
		tin_vec3 r2 = tin_sub_v3(c->position, b2->transform.translation);

		/* Relative velocity at contact */
		tin_vec3 v_rel = tin_sub_v3(
			tin_add_v3(b2->velocity, tin_cross_v3(b2->angular_velocity, r2)),
			tin_add_v3(b1->velocity, tin_cross_v3(b1->angular_velocity, r1)));

		/* Compute normal impulse */
		tin_scalar magnitude = c->normal_mass * (-tin_dot_v3(v_rel, c->normal) + c->bias);
		magnitude = MAX(magnitude, 0.0f);
		tin_vec3 impulse = tin_scale_v3(magnitude, c->normal);

		/* Apply contact impulse */
		b1->velocity = tin_sub_v3(b1->velocity, tin_scale_v3(b1->inv_mass, impulse));
		b1->angular_velocity = tin_sub_v3(b1->angular_velocity,
			tin_solve_inertia(b1, tin_cross_v3(r1, impulse)));

		b2->velocity = tin_add_v3(b2->velocity, tin_scale_v3(b2->inv_mass, impulse));
		b2->angular_velocity = tin_add_v3(b2->angular_velocity,
			tin_solve_inertia(b2, tin_cross_v3(r2, impulse)));
	}
}
