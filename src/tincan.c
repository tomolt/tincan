#include "tincan.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
	tin_vec3 former_dir = tin_apply_qt(tin_conjugate_qt(s->former_transform->rotation), dir);
	sup->rel_former = tin_polytope_support(s->former_polytope, former_dir);
	tin_vec3 former_abs = tin_apply_qt(s->former_transform->rotation, sup->rel_former);
	former_abs = tin_add_v3(s->former_transform->translation, former_abs);

	tin_vec3 latter_dir = tin_apply_qt(tin_conjugate_qt(s->latter_transform->rotation), tin_neg_v3(dir));
	sup->rel_latter = tin_polytope_support(s->latter_polytope, latter_dir);
	tin_vec3 latter_abs = tin_apply_qt(s->latter_transform->rotation, sup->rel_latter);
	latter_abs = tin_add_v3(s->latter_transform->translation, latter_abs);

	sup->abs = tin_sub_v3(former_abs, latter_abs);
}

int
tin_construct_portal(const tin_polysum *s, const tin_ray *r, tin_portal *p)
{
	tin_vec3 oa, ob, oc;

	tin_vec3 dir = r->dir;
	tin_polysum_support(s, dir, &p->a);
	oa = tin_sub_v3(p->a.abs, r->origin);

	dir = tin_sub_v3(r->dir, tin_normalize_v3(oa));
	if (tin_dot_v3(dir, dir) == 0.0f) {
		dir = tin_gram_schmidt(oa, (tin_vec3) {{ 1.0f, 0.0f, 0.0f }});
		dir = tin_normalize_v3(dir);
		if (tin_dot_v3(dir, dir) == 0.0f) {
			dir = (tin_vec3) {{ 0.0f, 0.0f, 1.0f }};
		}
	}

	tin_polysum_support(s, dir, &p->b);
	ob = tin_sub_v3(p->b.abs, r->origin);

	dir = tin_cross_v3(oa, ob);
	if (tin_dot_v3(dir, dir) == 0.0f) {
		return 0;
	}
	if (tin_dot_v3(dir, r->dir) < 0.0f) {
		dir = tin_neg_v3(dir);
		/* preserve ccw winding */
		tin_pspoint temp = p->a;
		p->a = p->b;
		p->b = temp;
	}
	
	for (int it = 0;; it++) {
		if (it >= 100) {
			fprintf(stderr, "Portal construction took too many iterations.\n");
			return 0;
		}

		if (tin_dot_v3(dir, dir) == 0.0f) {
			fprintf(stderr, "Portal collapsed during construction. (it=%d)\n", it);
			return 0;
		}

		tin_polysum_support(s, dir, &p->c);
		if (tin_dot_v3(dir, p->c.abs) <= 0.0f) {
			return 0;
		}
		
		tin_vec3 normal = tin_cross_v3(
			tin_sub_v3(p->b.abs, p->a.abs),
			tin_sub_v3(p->c.abs, p->a.abs));
		if (tin_dot_v3(r->dir, normal) < 0.0f) {
			fprintf(stderr, "Portal faces the wrong direction! (it=%d)\n", it);
		}

		oc = tin_sub_v3(p->c.abs, r->origin);
		tin_vec3 d_x_oc = tin_cross_v3(r->dir, oc);
		
		/* if <D , (OB x OC)> < 0 */
		ob = tin_sub_v3(p->b.abs, r->origin);
		if (-tin_dot_v3(ob, d_x_oc) < 0.0f) {
			p->a = p->c; /* throw away A */
			dir = tin_cross_v3(oc, ob);
			continue;
		}
		
		/* if <D , (OC x OA)> < 0 */
		oa = tin_sub_v3(p->a.abs, r->origin);
		if (tin_dot_v3(oa, d_x_oc) < 0.0f) {
			p->b = p->c; /* throw away B */
			dir = tin_cross_v3(oa, oc);
			continue;
		}
		
		return 1;
	}
}

#if 0

static int
refine_portal(const Minkowski *m, const Ray *r, Portal *p)
{
	for (int it = 0;; it++) {
		tin_vec3 ab = tin_sub_v3(p->b.global, p->a.global);
		tin_vec3 ac = tin_sub_v3(p->c.global, p->a.global);
		p->normal = tin_cross_v3(ab, ac);

		if (tin_dot_v3(p->normal, p->normal) == 0.0f) {
			log_warn("portal collapsed. (it=%d)", it);
			break;
		}

		if (it >= 100) {
			log_warn("MPR refine_portal() took too many iterations.");
			break;
		}

		SharedPoint s;
		minkowski_support(m, p->normal, &s);
		
		tin_vec3 as = tin_sub_v3(s.global, p->a.global);
		if (tin_dot_v3(as, p->normal) < 0.001f * tin_dot_v3(p->normal, p->normal)) {
			break;
		}
		
		tin_vec3 os = tin_sub_v3(s.global, r->origin);
		tin_vec3 d_x_os = tin_cross_v3(r->dir, os);
		
		/* if <D , (OS x OA)> > 0 */
		tin_vec3 oa = tin_sub_v3(p->a.global, r->origin);
		if (tin_dot_v3(oa, d_x_os) > 0.0f) {
			/* if <D , (OS x OB)> > 0 */
			tin_vec3 ob = tin_sub_v3(p->b.global, r->origin);
			if (tin_dot_v3(ob, d_x_os) > 0.0f) {
				p->a = s; /* SBC */
			} else {
				p->c = s; /* ABS */
			}
		} else {
			/* if <D , (OS x OC)> > 0 */
			tin_vec3 oc = tin_sub_v3(p->c.global, r->origin);
			if (tin_dot_v3(oc, d_x_os) > 0.0f) {
				p->b = s; /* ASC */
			} else {
				p->a = s; /* SBC */
			}
		}
	}

	if (tin_dot_v3(p->normal, p->a.global) >= 0.0f) {
		return 1;
	} else {
		return 0;
	}
}

static void
calculate_contact_point(const Minkowski *m, const Ray *r, const Portal *p, SharedPoint *point)
{
	tin_vec3 ab = tin_sub_v3(p->b.global, p->a.global);
	tin_vec3 ac = tin_sub_v3(p->c.global, p->a.global);
	float area_total = pgram_area(ab, ac);
	if (area_total == 0.0f) {
		*point = p->a;
		log_warn("Calculating contact point on collapsed portal.");
		return;
	}
	float proj = tin_dot_v3(p->normal, r->dir);
	if (proj == 0.0f) {
		*point = p->a;
		log_warn("Portal and contact normal have become parallel.");
		return;
	}
	float t = tin_dot_v3(p->normal, tin_sub_v3(p->a.global, r->origin)) / proj;
	tin_vec3 q = tin_saxpy_v3(t, r->dir, r->origin);
	tin_vec3 aq = tin_sub_v3(q, p->a.global);
	float area_b = pgram_area(aq, ac);
	float area_c = pgram_area(aq, ab);
	float beta = area_b / area_total;
	float gamma = area_c / area_total;
	float alpha = 1.0f - beta - gamma;
	
	if (alpha < 0.0f || alpha > 1.0001f
	||  beta  < 0.0f || beta  > 1.0001f
	||  gamma < 0.0f || gamma > 1.0001f) {
		log_warn("Barycentric coordinates lie outside of portal. (a=%f,b=%f,c=%f)", alpha, beta, gamma);
	}
	
	point->a_local = tin_scale_v3(alpha, p->a.a_local);
	point->a_local = tin_saxpy_v3(beta,  p->b.a_local, point->a_local);
	point->a_local = tin_saxpy_v3(gamma, p->c.a_local, point->a_local);

	point->b_local = tin_scale_v3(alpha, p->a.b_local);
	point->b_local = tin_saxpy_v3(beta,  p->b.b_local, point->b_local);
	point->b_local = tin_saxpy_v3(gamma, p->c.b_local, point->b_local);
}

int
mpr_intersect(const Minkowski *m, const Ray *ray, Contact *contact)
{
	Portal p;
	Ray r = *ray;
	if (!construct_portal(m, &r, &p)) {
		return 0;
	}
	if (!refine_portal(m, &r, &p)) {
		return 0;
	}
	p.normal = tin_normalize_v3(p.normal);

	for (int it = 0; it < 4; it++) {
		Ray nr;
		nr.dir    = p.normal;
		nr.origin = (tin_vec3) {{ 0.0f, 0.0f, 0.0f }};
		Portal np;
		if (!construct_portal(m, &nr, &np)) {
			break;
		}
		if (!refine_portal(m, &nr, &np)) {
			break;
		}
		np.normal = tin_normalize_v3(np.normal);

		float proj = tin_dot_v3(p.normal, np.normal);
		p = np; /* unintentional, i swear */
		r = nr;
		if (proj >= 0.99f) break;
	}

	SharedPoint point;
	calculate_contact_point(m, &r, &p, &point);
	contact->normal   = p.normal;
	contact->position = point.global;
	return 1;
}

#endif
