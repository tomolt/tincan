#include "tincan.h"

#include <string.h>
#include <stdio.h>
#include <math.h>

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

Tin_Vec3
tin_neg_v3(Tin_Vec3 x)
{
	return (Tin_Vec3) {{ -x.c[0], -x.c[1], -x.c[2] }};
}

Tin_Vec3
tin_add_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return (Tin_Vec3) {{ a.c[0]+b.c[0], a.c[1]+b.c[1], a.c[2]+b.c[2] }};
}

Tin_Vec3
tin_sub_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return (Tin_Vec3) {{ a.c[0]-b.c[0], a.c[1]-b.c[1], a.c[2]-b.c[2] }};
}

Tin_Vec3
tin_scale_v3(Tin_Scalar a, Tin_Vec3 x)
{
	return (Tin_Vec3) {{ a*x.c[0], a*x.c[1], a*x.c[2] }};
}

Tin_Vec3
tin_saxpy_v3(Tin_Scalar a, Tin_Vec3 x, Tin_Vec3 y)
{
	return tin_add_v3(tin_scale_v3(a, x), y);
}

Tin_Vec3
tin_cross_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	Tin_Vec3 c;
	c.c[0] = a.c[1] * b.c[2] - a.c[2] * b.c[1];
	c.c[1] = a.c[2] * b.c[0] - a.c[0] * b.c[2];
	c.c[2] = a.c[0] * b.c[1] - a.c[1] * b.c[0];
	return c;
}

Tin_Scalar
tin_dot_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return a.c[0]*b.c[0] + a.c[1]*b.c[1] + a.c[2]*b.c[2];
}

Tin_Scalar
tin_length_v3(Tin_Vec3 v)
{
	return sqrtf(tin_dot_v3(v, v));
}

Tin_Vec3
tin_normalize_v3(Tin_Vec3 v)
{
	Tin_Scalar norm = tin_length_v3(v);
	if (norm > 0.0f) {
		return tin_scale_v3(1.0f / norm, v);
	} else {
		return v;
	}
}

Tin_Vec3
tin_hadamard_v3(Tin_Vec3 a, Tin_Vec3 b)
{
	return (Tin_Vec3) {{ a.c[0]*b.c[0], a.c[1]*b.c[1], a.c[2]*b.c[2] }};
}

Tin_Quat
tin_make_qt(Tin_Vec3 axis, Tin_Scalar angle)
{
	Tin_Scalar h = angle / 2.0f;
	return (Tin_Quat) { tin_scale_v3(sinf(h), axis), cosf(h) };
}

Tin_Quat
tin_mul_qt(Tin_Quat a, Tin_Quat b)
{
	Tin_Quat d;
	d.v = tin_cross_v3(a.v, b.v);
	d.v = tin_saxpy_v3(a.s, b.v, d.v);
	d.v = tin_saxpy_v3(b.s, a.v, d.v);
	d.s = a.s * b.s - tin_dot_v3(a.v, b.v);
	return d;
}

Tin_Vec3
tin_apply_qt(Tin_Quat q, Tin_Vec3 v)
{
	Tin_Vec3 t = tin_cross_v3(q.v, v);
	Tin_Vec3 d = tin_cross_v3(q.v, t);
	d = tin_saxpy_v3(q.s, t, d);
	d = tin_saxpy_v3(2.0f, d, v);
	return d;
}

Tin_Quat
tin_conjugate_qt(Tin_Quat q)
{
	return (Tin_Quat) { tin_neg_v3(q.v), q.s };
}

Tin_Vec3
tin_gram_schmidt(Tin_Vec3 fixed, Tin_Vec3 var)
{
	Tin_Scalar fscale = tin_dot_v3(fixed, fixed);
	if (fscale) {
		return tin_saxpy_v3(-tin_dot_v3(fixed, var) / fscale, fixed, var);
	} else {
		return var;
	}
}

Tin_Scalar
tin_prlgram_area(Tin_Vec3 e1, Tin_Vec3 e2)
{
	Tin_Vec3 perp = tin_gram_schmidt(e1, e2);
	return sqrtf(tin_dot_v3(e1, e1) * tin_dot_v3(perp, perp));
}

Tin_Vec3
tin_fwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec)
{
	vec = tin_scale_v3(transform->scale,       vec);
	vec = tin_apply_qt(transform->rotation,    vec);
	vec = tin_add_v3  (transform->translation, vec);
	return vec;
}

Tin_Vec3
tin_fwtrf_dir(const Tin_Transform *transform, Tin_Vec3 vec)
{
	return tin_apply_qt(transform->rotation, vec);
}

Tin_Vec3
tin_bwtrf_point(const Tin_Transform *transform, Tin_Vec3 vec)
{
	vec = tin_sub_v3  (vec, transform->translation);
	vec = tin_apply_qt(tin_conjugate_qt(transform->rotation), vec);
	vec = tin_scale_v3(1.0f / transform->scale, vec);
	return vec;
}

Tin_Vec3
tin_bwtrf_dir(const Tin_Transform *transform, Tin_Vec3 vec)
{
	return tin_apply_qt(tin_conjugate_qt(transform->rotation), vec);
}

Tin_Vec3
tin_polytope_support(const Tin_Polytope *polytope, Tin_Vec3 dir)
{
	Tin_Scalar bestScore = -INFINITY;
	int bestIdx = -1;
	for (int idx = 0; idx < polytope->numVertices; idx++) {
		Tin_Scalar score = tin_dot_v3(polytope->vertices[idx], dir);
		if (score > bestScore) {
			bestScore = score;
			bestIdx   = idx;
		}
	}
	return polytope->vertices[bestIdx];
}

void
tin_polysum_support(const Tin_Polysum *s, Tin_Vec3 dir, Tin_Pspoint *sup)
{
	Tin_Vec3 former_dir = tin_bwtrf_dir(s->transform1, dir);
	sup->relTo1 = tin_polytope_support(s->polytope1, former_dir);
	Tin_Vec3 former_abs = tin_fwtrf_point(s->transform1, sup->relTo1);

	Tin_Vec3 latter_dir = tin_bwtrf_dir(s->transform2, tin_neg_v3(dir));
	sup->relTo2 = tin_polytope_support(s->polytope2, latter_dir);
	Tin_Vec3 latter_abs = tin_fwtrf_point(s->transform2, sup->relTo2);

	sup->abs = tin_sub_v3(former_abs, latter_abs);
}

int
tin_construct_portal(const Tin_Polysum *ps, const Tin_Ray *r, Tin_Portal *p)
{
	/* find the very first point */
	tin_polysum_support(ps, r->dir, &p->a);

	/* find the second point */
	{
		Tin_Vec3 oa = tin_sub_v3(p->a.abs, r->origin);
		Tin_Vec3 dir = tin_gram_schmidt(oa, r->dir);
		if (tin_dot_v3(dir, dir) == 0.0f) {
			dir = tin_gram_schmidt(oa, (Tin_Vec3) {{ 1.0f, 0.0f, 0.0f }});
			if (tin_dot_v3(dir, dir) == 0.0f) {
				dir = (Tin_Vec3) {{ 0.0f, 0.0f, 1.0f }};
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
		Tin_Vec3 oa = tin_sub_v3(p->a.abs, r->origin);
		Tin_Vec3 ob = tin_sub_v3(p->b.abs, r->origin);
		Tin_Vec3 dir = tin_cross_v3(oa, ob);
		if (tin_dot_v3(dir, r->dir) < 0.0f) {
			dir = tin_neg_v3(dir);
		}
		tin_polysum_support(ps, dir, &p->c);
		Tin_Scalar score = tin_dot_v3(dir, p->c.abs);
		if (score <= tin_dot_v3(dir, p->a.abs) || score <= tin_dot_v3(dir, p->b.abs)) {
			return 0;
		}
		Tin_Vec3 oc = tin_sub_v3(p->c.abs, r->origin);
		
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
tin_refine_portal(const Tin_Polysum *ps, const Tin_Ray *r, Tin_Portal *p)
{
	Tin_Vec3 ab = tin_sub_v3(p->b.abs, p->a.abs);
	Tin_Vec3 ac = tin_sub_v3(p->c.abs, p->a.abs);
	p->normal = tin_cross_v3(ab, ac);
	if (tin_dot_v3(p->normal, r->dir) < 0.0f) {
		Tin_Pspoint temp;
		temp = p->a;
		p->a = p->b;
		p->b = temp;
	}

	for (int it = 0;; it++) {
		Tin_Vec3 ab = tin_sub_v3(p->b.abs, p->a.abs);
		Tin_Vec3 ac = tin_sub_v3(p->c.abs, p->a.abs);
		p->normal = tin_cross_v3(ab, ac);

		if (tin_dot_v3(p->normal, p->normal) == 0.0f) {
			fprintf(stderr, "portal collapsed. (it=%d)\n", it);
			break;
		}

		if (it >= 100) {
			fprintf(stderr, "MPR refine_portal() took too many iterations.\n");
			break;
		}

		Tin_Pspoint s;
		tin_polysum_support(ps, p->normal, &s);
		
		{
			Tin_Vec3 as = tin_sub_v3(s.abs, p->a.abs);
			Tin_Vec3 bs = tin_sub_v3(s.abs, p->b.abs);
			Tin_Vec3 cs = tin_sub_v3(s.abs, p->c.abs);
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
		
		Tin_Vec3 os = tin_sub_v3(s.abs, r->origin);
		Tin_Vec3 d_x_os = tin_cross_v3(r->dir, os);
		
		/* if <D , (OS x OA)> > 0 */
		Tin_Vec3 oa = tin_sub_v3(p->a.abs, r->origin);
		if (tin_dot_v3(oa, d_x_os) > 0.0f) {
			/* if <D , (OS x OB)> > 0 */
			Tin_Vec3 ob = tin_sub_v3(p->b.abs, r->origin);
			if (tin_dot_v3(ob, d_x_os) > 0.0f) {
				p->a = s; /* SBC */
			} else {
				p->c = s; /* ABS */
			}
		} else {
			/* if <D , (OS x OC)> > 0 */
			Tin_Vec3 oc = tin_sub_v3(p->c.abs, r->origin);
			if (tin_dot_v3(oc, d_x_os) > 0.0f) {
				p->b = s; /* ASC */
			} else {
				p->a = s; /* SBC */
			}
		}
	}
}

void
tin_calculate_contact(const Tin_Ray *r, const Tin_Portal *p, Tin_Contact *contact)
{
	contact->normal = tin_normalize_v3(p->normal);

	Tin_Vec3 ab = tin_sub_v3(p->b.abs, p->a.abs);
	Tin_Vec3 ac = tin_sub_v3(p->c.abs, p->a.abs);
	Tin_Scalar area_total = tin_prlgram_area(ab, ac);
	if (area_total == 0.0f) {
		contact->rel1 = p->a.relTo1;
		contact->rel2 = p->a.relTo2;
		fprintf(stderr, "Calculating contact point on collapsed portal.\n");
		return;
	}
	Tin_Scalar proj = tin_dot_v3(p->normal, r->dir);
	if (proj == 0.0f) {
		contact->rel1 = p->a.relTo1;
		contact->rel2 = p->a.relTo2;
		fprintf(stderr, "Portal and contact normal have become parallel.\n");
		return;
	}
	Tin_Scalar t = tin_dot_v3(p->normal, tin_sub_v3(p->a.abs, r->origin)) / proj;
	Tin_Vec3 q = tin_saxpy_v3(t, r->dir, r->origin);
	Tin_Vec3 aq = tin_sub_v3(q, p->a.abs);
	Tin_Scalar area_b = tin_prlgram_area(aq, ac);
	Tin_Scalar area_c = tin_prlgram_area(aq, ab);
	Tin_Scalar beta = area_b / area_total;
	Tin_Scalar gamma = area_c / area_total;
	Tin_Scalar alpha = 1.0f - beta - gamma;
	
	contact->rel1 = tin_scale_v3(alpha, p->a.relTo1);
	contact->rel1 = tin_saxpy_v3(beta,  p->b.relTo1, contact->rel1);
	contact->rel1 = tin_saxpy_v3(gamma, p->c.relTo1, contact->rel1);

	contact->rel2 = tin_scale_v3(alpha, p->a.relTo2);
	contact->rel2 = tin_saxpy_v3(beta,  p->b.relTo2, contact->rel2);
	contact->rel2 = tin_saxpy_v3(gamma, p->c.relTo2, contact->rel2);
}

int
tin_polytope_collide(
	const Tin_Polytope *pa, const Tin_Transform *ta,
	const Tin_Polytope *pb, const Tin_Transform *tb,
	Tin_Contact *contact)
{
	Tin_Polysum ps = { pa, ta, pb, tb };

	Tin_Ray r;
	r.origin = tin_sub_v3(ta->translation, tb->translation);
	r.dir    = tin_neg_v3(r.origin);
	Tin_Scalar norm = sqrtf(tin_dot_v3(r.dir, r.dir));
	if (norm == 0.0f) {
		/* FIXME */
		return 1;
	}
	r.dir = tin_normalize_v3(r.dir);

	Tin_Portal p;
	if (!tin_construct_portal(&ps, &r, &p)) {
		return 0;
	}
	tin_refine_portal(&ps, &r, &p);
	p.normal = tin_normalize_v3(p.normal);
	if (tin_dot_v3(p.normal, p.a.abs) <= 0.0f) {
		return 0;
	}

	for (int it = 0; it < 4; it++) {
		Tin_Ray nr;
		nr.dir    = p.normal;
		nr.origin = (Tin_Vec3) {{ 0.0f, 0.0f, 0.0f }};
		Tin_Portal np;
		if (!tin_construct_portal(&ps, &nr, &np)) {
			break;
		}
		tin_refine_portal(&ps, &nr, &np);
		np.normal = tin_normalize_v3(np.normal);

		Tin_Scalar proj = tin_dot_v3(p.normal, np.normal);
		p = np; /* unintentional, i swear */
		r = nr;
		if (proj >= 0.99f) break;
	}

	tin_calculate_contact(&r, &p, contact);
	return 1;
}

Tin_Vec3
tin_solve_inertia(const Tin_Body *body, Tin_Vec3 vec)
{
	vec = tin_bwtrf_dir(&body->transform, vec);
	vec = tin_hadamard_v3(body->invInertia, vec);
	vec = tin_fwtrf_dir(&body->transform, vec);
	return vec;
}

void
tin_arbiter_update(Tin_Arbiter *arbiter)
{
	const Tin_Scalar maxSeparation = 0.1f;
	const Tin_Scalar maxStretch = 0.3f;

	for (int i = 0; i < arbiter->numContacts; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		Tin_Vec3 p1 = tin_fwtrf_point(&arbiter->body1->transform, contact->rel1);
		Tin_Vec3 p2 = tin_fwtrf_point(&arbiter->body2->transform, contact->rel2);

		Tin_Scalar separation = tin_dot_v3(contact->normal, tin_sub_v3(p2, p1));

		if (separation > maxSeparation) {
			*contact = arbiter->contacts[--arbiter->numContacts];
			i--;
			continue;
		}

		Tin_Scalar stretch = tin_length_v3(tin_gram_schmidt(contact->normal, tin_sub_v3(p1, p2)));
		if (stretch > maxStretch) {
			*contact = arbiter->contacts[--arbiter->numContacts];
			i--;
			continue;
		}
	}

}

void
tin_arbiter_add_contact(Tin_Arbiter *arbiter, Tin_Contact contact)
{
	Tin_Vec3 p1, p2;
	p1 = tin_fwtrf_point(&arbiter->body1->transform, contact.rel1);
	p2 = tin_fwtrf_point(&arbiter->body2->transform, contact.rel2);
	contact.baseStretch = tin_length_v3(tin_gram_schmidt(contact.normal, tin_sub_v3(p1, p2)));

	if (arbiter->numContacts < TIN_MAX_CONTACTS) {
		arbiter->contacts[arbiter->numContacts++] = contact;
		return;
	}

	Tin_Vec3 newPos = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
	int bestIdx = -1;
	Tin_Scalar bestDist = INFINITY;
	for (int idx = 0; idx < arbiter->numContacts; idx++) {
		Tin_Contact *c = &arbiter->contacts[idx];
		p1 = tin_fwtrf_point(&arbiter->body1->transform, c->rel1);
		p2 = tin_fwtrf_point(&arbiter->body2->transform, c->rel2);
		Tin_Vec3 oldPos = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
		Tin_Vec3 diff = tin_sub_v3(oldPos, newPos);
		Tin_Scalar dist = tin_dot_v3(diff, diff);
		if (dist < bestDist) {
			bestIdx  = idx;
			bestDist = dist;
		}
	}
	arbiter->contacts[bestIdx] = contact;
}

void
tin_arbiter_prestep(Tin_Arbiter *arbiter, Tin_Scalar invDt)
{
	const Tin_Scalar allowedPenetration = 0.01f;
	const Tin_Scalar biasFactor = 0.1f;

	for (int i = 0; i < arbiter->numContacts; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];

		Tin_Vec3 p1 = tin_fwtrf_point(&arbiter->body1->transform, contact->rel1);
		Tin_Vec3 p2 = tin_fwtrf_point(&arbiter->body2->transform, contact->rel2);

		contact->position   = tin_scale_v3(0.5f, tin_add_v3(p1, p2));
		contact->separation = tin_dot_v3(contact->normal, tin_sub_v3(p2, p1));

		Tin_Vec3 r1 = tin_sub_v3(contact->position, arbiter->body1->transform.translation);
		Tin_Vec3 r2 = tin_sub_v3(contact->position, arbiter->body2->transform.translation);

		/* Precompute normal mass and bias */
		contact->normalMass = 1.0f / (arbiter->body1->invMass + arbiter->body2->invMass + tin_dot_v3(contact->normal, tin_add_v3(
			tin_cross_v3(tin_solve_inertia(arbiter->body1, tin_cross_v3(r1, contact->normal)), r1),
			tin_cross_v3(tin_solve_inertia(arbiter->body2, tin_cross_v3(r2, contact->normal)), r2))));

		contact->bias = -biasFactor * invDt * MIN(0.0f, contact->separation + allowedPenetration);
	}
}

void
tin_apply_impulse(Tin_Body *body, Tin_Vec3 impulse, Tin_Vec3 at)
{
	body->velocity = tin_add_v3(body->velocity, tin_scale_v3(body->invMass, impulse));
	body->angular_velocity = tin_add_v3(body->angular_velocity,
		tin_solve_inertia(body, tin_cross_v3(at, impulse)));
}

void
tin_arbiter_apply_impulse(Tin_Arbiter *arbiter, Tin_Scalar invDt)
{
	for (int i = 0; i < arbiter->numContacts; i++) {
		Tin_Contact *contact = &arbiter->contacts[i];
		if (contact->separation > 0.0f) continue;

		Tin_Vec3 r1 = tin_sub_v3(contact->position, arbiter->body1->transform.translation);
		Tin_Vec3 r2 = tin_sub_v3(contact->position, arbiter->body2->transform.translation);

		/* Relative velocity at contact */
		Tin_Vec3 v_rel = tin_sub_v3(
			tin_add_v3(arbiter->body2->velocity, tin_cross_v3(arbiter->body2->angular_velocity, r2)),
			tin_add_v3(arbiter->body1->velocity, tin_cross_v3(arbiter->body1->angular_velocity, r1)));

		/* Compute normal impulse */
		Tin_Scalar magnitude = contact->normalMass * (-tin_dot_v3(v_rel, contact->normal) + contact->bias);
		magnitude = MAX(magnitude, 0.0f);
		Tin_Vec3 impulse = tin_scale_v3(magnitude, contact->normal);

		/* Apply contact impulse */
		tin_apply_impulse(arbiter->body1, tin_neg_v3(impulse), r1);
		tin_apply_impulse(arbiter->body2, impulse, r2);

		/* Compute friction impulse */
		Tin_Scalar friction_coefficient = 0.2f;
		Tin_Vec3 friction = tin_gram_schmidt(contact->normal, v_rel);
		if (tin_dot_v3(friction, friction) <= 10.0f * TIN_EPSILON) {
			friction_coefficient *= 2.0f;
		}
		friction = tin_scale_v3(-friction_coefficient / invDt, friction);
		if (arbiter->body1->invMass > 0.0f) {
			tin_apply_impulse(arbiter->body1, tin_scale_v3(-1.0f / arbiter->body1->invMass, friction), r1);
		}
		if (arbiter->body2->invMass > 0.0f) {
			tin_apply_impulse(arbiter->body2, tin_scale_v3(1.0f / arbiter->body2->invMass, friction), r2);
		}
	}
}

Tin_Arbiter *
tin_find_arbiter(Tin_Scene *scene, Tin_Body *body1, Tin_Body *body2)
{
	TIN_FOR_EACH(arbiter, scene->arbiters, Tin_Arbiter, node) {
		if (arbiter->body1 == body1 && arbiter->body2 == body2) {
			return arbiter;
		}
	}
	arbiter = tin_add_arbiter(scene);
	arbiter->body1 = body1;
	arbiter->body2 = body2;
	return arbiter;
}

void
tin_scene_update(Tin_Scene *scene)
{
	TIN_FOR_EACH(arbiter, scene->arbiters, Tin_Arbiter, node) {
		tin_arbiter_update(arbiter);
	}
}

void
tin_check_collision(Tin_Scene *scene, Tin_Body *body1, Tin_Body *body2)
{
	Tin_Arbiter *arbiter = tin_find_arbiter(scene, body1, body2);
	if (!arbiter) return;

	Tin_Contact contact;
	int colliding = tin_polytope_collide(
			body1->shapeParams, &body1->transform,
			body2->shapeParams, &body2->transform, &contact);

	if (colliding) {
		printf("D %f\n", tin_dot_v3(contact.normal,
					tin_sub_v3(body2->transform.translation, body1->transform.translation)));
		tin_arbiter_add_contact(arbiter, contact);
	}
}

void
tin_scene_prestep(Tin_Scene *scene, Tin_Scalar invDt)
{
	TIN_FOR_EACH(arbiter, scene->arbiters, Tin_Arbiter, node) {
		tin_arbiter_prestep(arbiter, invDt);
	}

}

void
tin_scene_step(Tin_Scene *scene, Tin_Scalar invDt)
{
	TIN_FOR_EACH(arbiter, scene->arbiters, Tin_Arbiter, node) {
		tin_arbiter_apply_impulse(arbiter, invDt);
	}
}

void
tin_integrate(Tin_Scene *scene, Tin_Scalar dt)
{
	TIN_FOR_EACH(body, scene->bodies, Tin_Body, node) {
		if (tin_dot_v3(body->velocity, body->velocity) > 10000.0f * TIN_EPSILON) {
			body->transform.translation = tin_saxpy_v3(dt, body->velocity, body->transform.translation);
		}
		
		Tin_Vec3 av = body->angular_velocity;
		Tin_Scalar angle = sqrtf(tin_dot_v3(av, av));
		
		if (fabs(angle) > 100.0f * TIN_EPSILON) {
			body->transform.rotation = tin_mul_qt(
				tin_make_qt(tin_normalize_v3(av), angle * dt), body->transform.rotation);
		}

		if (body->invMass != 0.0f) {
			body->velocity = tin_saxpy_v3(dt, (Tin_Vec3) {{ 0.0f, -6.0f, 0.0f }}, body->velocity);
		}
	}
}

void
tin_broadphase(Tin_Scene *scene, void (*func)(Tin_Scene *, Tin_Body *, Tin_Body *))
{
	TIN_FOR_EACH(body1, scene->bodies, Tin_Body, node) {
		const Tin_Polytope *polytope1 = body1->shapeParams;
		TIN_FOR_RANGE(body2, body1->node, scene->bodies, Tin_Body, node) {
			const Tin_Polytope *polytope2 = body2->shapeParams;
			Tin_Vec3 diff = tin_sub_v3(body2->transform.translation, body1->transform.translation);
			Tin_Scalar radii =
				polytope2->radius * body2->transform.scale +
				polytope1->radius * body1->transform.scale;
			if (tin_dot_v3(diff, diff) <= radii * radii) {
				func(scene, body1, body2);
			}
		}
	}
}

void
tin_simulate(Tin_Scene *scene, Tin_Scalar dt)
{
	Tin_Scalar stepDt = dt / 4.0f;
	Tin_Scalar stepInvDt = 1.0f / stepDt;
	for (int step = 0; step < 4; step++) {
		tin_scene_update(scene);
		tin_broadphase(scene, tin_check_collision);
		tin_scene_prestep(scene, stepInvDt);
		for (int iter = 0; iter < 4; iter++) {
			tin_scene_step(scene, stepInvDt);
		}
		tin_integrate(scene, stepDt);
	}
}

Tin_Body *
tin_add_body(Tin_Scene *scene)
{
	Tin_Body *body = scene->bodyAllocator.alloc(scene->bodyAllocator.userPointer);
	memset(body, 0, sizeof *body);
	TIN_LIST_LINK(*scene->bodies.prev, body->node);
	TIN_LIST_LINK(body->node, scene->bodies);
	return body;
}

Tin_Arbiter *
tin_add_arbiter(Tin_Scene *scene)
{
	Tin_Arbiter *arbiter = scene->arbiterAllocator.alloc(scene->arbiterAllocator.userPointer);
	memset(arbiter, 0, sizeof *arbiter);
	TIN_LIST_LINK(*scene->arbiters.prev, arbiter->node);
	TIN_LIST_LINK(arbiter->node, scene->arbiters);
	return arbiter;
}

