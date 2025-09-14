#include <tincan.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include "mat4.h"
#include "render.h"
#include "camera.h"

#define MAX_OBJECTS 100

typedef struct {
	Tin_Body *body;
	Model    *model;
	Tin_Vec3  color;
} Object;

Tin_Scene scene;
Tin_SweepPrune sweep_prune;

static Camera camera;
static Object objects[MAX_OBJECTS];
static int    num_objects;

static Tin_Shape cone_shape;
static Model     cone_model;
static Tin_Shape cube_shape;
static Model     cube_model;

static float time_multiplier = 1.0f;

static bool throwing;
static double latest_throw;

static void
init_cone(int tessel)
{
	int nverts   = tessel + 2;
	int nindices = 6 * tessel;
	Tin_Vec3 *verts   = calloc(nverts,   sizeof (Tin_Vec3));
	GLushort *indices = calloc(nindices, sizeof (GLushort));

	for (int v = 0; v < tessel; v++) {
		float angle = 2.0f * M_PI * v / tessel;
		verts[v] = (Tin_Vec3){{ cosf(angle), 0.5f, sinf(angle) }};
	}
	verts[tessel + 0] = (Tin_Vec3){{ 0.0f, -1.5f, 0.0f }};
	verts[tessel + 1] = (Tin_Vec3){{ 0.0f, 0.5f, 0.0f }};

	for (int v = 0; v < tessel; v++) {
		int w = (v + 1) % tessel;
		indices[6*v+0] = v;
		indices[6*v+1] = w;
		indices[6*v+2] = tessel + 0;
		indices[6*v+3] = w;
		indices[6*v+4] = v;
		indices[6*v+5] = tessel + 1;
	}

	cone_shape.kind = TIN_POLYTOPE;
	cone_shape.polytope.vertices = verts;
	cone_shape.polytope.numVertices = nverts;
	cone_shape.radius = 1.5f;
	Tin_Scalar r = 1.0f, h = 2.0f;
	cone_shape.invInertia = (Tin_Vec3) {{
		1.0 / (3.0 / 80.0 * (h * h + 4.0 * r * r)),
		1.0 / (3.0 / 10.0),
		1.0 / (3.0 / 80.0 * (h * h + 4.0 * r * r)),
	}};

	cone_model = render_make_model(nverts, verts, nindices, indices);
	
	free(indices);
}

static void
init_cube(void)
{
	static Tin_Vec3 verts[] = {
		{{ -1, -1, -1 }},
		{{  1, -1, -1 }},
		{{ -1,  1, -1 }},
		{{  1,  1, -1 }},
		{{ -1, -1,  1 }},
		{{  1, -1,  1 }},
		{{ -1,  1,  1 }},
		{{  1,  1,  1 }},
	};
	static GLushort indices[] = {
		0, 2, 3, 3, 1, 0,
		4, 5, 7, 7, 6, 4,
		0, 1, 5, 5, 4, 0,
		2, 6, 7, 7, 3, 2,
		0, 4, 6, 6, 2, 0,
		1, 3, 7, 7, 5, 1,
	};
	static int faceIndices[] = {
		0, 2, 3, 1,
		4, 5, 7, 6,
		0, 1, 5, 4,
		2, 6, 7, 3,
		0, 4, 6, 2,
		1, 3, 7, 5,
	};
	static int faceOffsets[] = {
		0,
		4,
		8,
		12,
		16,
		20,
		24,
	};
	static Tin_Vec3 faceNormals[] = {
		{{ 0.0, 0.0,-1.0}},
		{{ 0.0, 0.0, 1.0}},
		{{ 0.0,-1.0, 0.0}},
		{{ 0.0, 1.0, 0.0}},
		{{-1.0, 0.0, 0.0}},
		{{ 1.0, 0.0, 0.0}},
	};

	cube_shape.kind = TIN_POLYTOPE;
	cube_shape.polytope.vertices = verts;
	cube_shape.polytope.numVertices = 8;
	cube_shape.polytope.faceIndices = faceIndices;
	cube_shape.polytope.faceOffsets = faceOffsets;
	cube_shape.polytope.faceNormals = faceNormals;
	cube_shape.polytope.numFaces = 6;
	cube_shape.radius = sqrtf(3.0f);
	cube_shape.invInertia = (Tin_Vec3) {{ 6.0f, 6.0f, 6.0f }};
	
	cube_model = render_make_model(8, verts, 6 * 6, indices);
}

static Camera *cam = &camera;

static void
key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	(void)scancode, (void)mods;

	if (action == GLFW_REPEAT) return;
	
	switch (key) {
	case GLFW_KEY_ESCAPE:
		if (action == GLFW_PRESS) {
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
		break;

	case GLFW_KEY_W: cam->f = action == GLFW_PRESS; break;
	case GLFW_KEY_A: cam->l = action == GLFW_PRESS; break;
	case GLFW_KEY_S: cam->b = action == GLFW_PRESS; break;
	case GLFW_KEY_D: cam->r = action == GLFW_PRESS; break;
	
	case GLFW_KEY_E:
		if (action == GLFW_RELEASE) {
			int mode = glfwGetInputMode(window, GLFW_CURSOR);
			if (mode == GLFW_CURSOR_NORMAL) {
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
				glfwGetCursorPos(window, &cam->cursor_x, &cam->cursor_y);
			} else {
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			}
		}
		break;

	case GLFW_KEY_B:
		time_multiplier = 1.0f;
		break;
	
	case GLFW_KEY_M:
		time_multiplier *= 2.0f;
		break;
	case GLFW_KEY_N:
		time_multiplier /= 2.0f;
		break;

	case GLFW_KEY_SPACE:
		if (action == GLFW_PRESS) {
			throwing = true;
			latest_throw = -INFINITY;
		} else {
			throwing = false;
		}
		break;
	}
}

static void
mouse_callback(GLFWwindow *window, double x, double y)
{
	(void)window;

	int mode = glfwGetInputMode(window, GLFW_CURSOR);
	if (mode != GLFW_CURSOR_DISABLED) return;

        float dx = x - cam->cursor_x;
        float dy = y - cam->cursor_y;

        cam->yaw -= dx * 0.001f;
        while (cam->yaw < 0.0f) {
                cam->yaw += 2.0f * M_PI;
        }
        while (cam->yaw >= 2.0f * M_PI) {
                cam->yaw -= 2.0f * M_PI;
        }

        cam->pitch -= dy * 0.001f;
        if (cam->pitch < -0.5f * M_PI) {
                cam->pitch = -0.5f * M_PI;
        }
        if (cam->pitch > 0.5f * M_PI) {
                cam->pitch = 0.5f * M_PI;
        }

        cam->cursor_x = x;
        cam->cursor_y = y;
}

static void
draw_objects(void)
{
	for (int o = 0; o < num_objects; o++) {
		const Object *obj = &objects[o];
		Tin_Vec3 color = obj->color;
		if (tin_island_find(obj->body)->islandStable && obj->body->invMass != 0.0) {
			color = tin_saxpy_v3(0.5, color, TIN_VEC3(0.5, 0.5, 0.5));
		}
		render_draw_model(obj->model, &obj->body->transform, color);
	}
}

void *
custom_alloc(void *userPointer)
{
	return malloc((size_t)userPointer);
}

void
custom_free(void *userPointer, void *memoryPointer)
{
	(void)userPointer;
	free(memoryPointer);
}

int
main(void)
{
	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow *window = glfwCreateWindow(800, 600, "tincan physics demo", NULL, NULL);
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, mouse_callback);

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	gladLoadGL(glfwGetProcAddress);
	glEnable(GL_DEPTH_TEST);
	render_init();
	render_light_dir = tin_normalize_v3(TIN_VEC3(-0.3, -0.6, 0.0));
	init_cone(16);
	init_cube();

	memcpy(camera.rotation, (Tin_Scalar[]){ 1, 0, 0, 0, 1, 0, 0, 0, 1 }, sizeof camera.rotation);
	glfwGetCursorPos(window, &camera.cursor_x, &camera.cursor_y);
	camera.position.c[2] = 5.0f;

	tin_create_pairtable(&scene.contactCache);
	tin_create_sweep_prune(&sweep_prune, &scene);
	scene.sweepPrune = &sweep_prune;

	/*
	Tin_Body *body1 = tin_add_body(&scene, &cone_shape, 1.0 / 1.0);
	body1->transform.translation = TIN_VEC3(0.0, 1.5, 0.0);
	objects[num_objects++] = (Object){
		body1,
		&cone_model,
		{{ 0.5f, 1.0f, 0.5f }}
	};
	*/

	Tin_Body *body2 = tin_add_body(&scene, &cube_shape, 1.0 / 3.0);
	body2->transform.translation = TIN_VEC3(0.0, -1.0, 0.0);
	body2->velocity = TIN_VEC3(-3, 1, 0);
	objects[num_objects++] = (Object){
		body2,
		&cube_model,
		{{ 1.0f, 0.5f, 0.5f }}
	};

	Tin_Body *body3 = tin_add_body(&scene, &cube_shape, 0.0);
	body3->transform.translation = TIN_VEC3(0.0, -23.0, 0.0);
	body3->transform.scale = 20.0;
	objects[num_objects++] = (Object){
		body3,
		&cube_model,
		{{ 1.0f, 1.0f, 1.0f }}
	};

	double accumTimings[6];

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		if (throwing) {
			if (glfwGetTime() - latest_throw > 0.2) {
				latest_throw = glfwGetTime();
				if (num_objects >= MAX_OBJECTS) {
					printf("Cannot spawn box: too many objects already present.\n");
				} else {
					Tin_Vec3 forward = {{ -cam->rotation[6], -cam->rotation[7], -cam->rotation[8] }};
					Tin_Body *body = tin_add_body(&scene, &cube_shape, 1.0 / 1.0);
					memcpy(body->transform.rotation, cam->rotation, sizeof cam->rotation);
					body->transform.translation = tin_saxpy_v3(1.0, forward, cam->position);
					body->transform.scale = 0.2;
					body->velocity = tin_scale_v3(5.0, forward);
					objects[num_objects++] = (Object) {
						body,
						&cube_model,
						{{ 0.5f, 0.5f, 1.0f }}
					};
				}
			}
		}

		Tin_Scalar dt = 1.0f / 60.0f * time_multiplier;

		camera_update(&camera, 1.0f / 60.0f);

		int width, height;
		glfwGetWindowSize(window, &width, &height);

		Tin_Transform camtrf;
		camera_transform(&camera, &camtrf);

		Tin_Transform ident = { .rotation = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }, .scale = 1.0f };

		double startTime = glfwGetTime();
		double timings[6];
		tin_simulate(&scene, dt, glfwGetTime, timings);
		double elapsedTime = glfwGetTime() - startTime;
		for (int t = 0; t < 6; t++) {
			accumTimings[t] *= 0.8;
			accumTimings[t] += 0.2 * timings[t];
		}

		char title[100];
		snprintf(title, sizeof title, "tincan physics demo - %02dms", (int)(elapsedTime * 1000.0f + 0.5f));
		glfwSetWindowTitle(window, title);

		render_proj_matrix = mat4_orthographic(-25.0, 25.0, -25.0, 25.0, 1.0, 50.0);
		render_view_matrix = mat4_look_at(tin_scale_v3(-25.0, render_light_dir), TIN_VEC3(0, 0, 0), TIN_VEC3(0.0, 1.0, 0.0));
		render_light_matrix = mat4_multiply(render_proj_matrix, render_view_matrix);
		render_start_shadow_pass();
		draw_objects();

		render_proj_matrix = mat4_perspective(70.0f, (float) width / height, 1.0f, 100.0f);
		render_view_matrix = mat4_from_inverse_transform(&camtrf);
		render_start_scene_pass(width, height);
		draw_objects();

#if 0
		render_start_overlay();
		for (int o = 0; o < num_objects; o++) {
			Object *obj = &objects[o];
			Tin_Body *body = obj->body;
			Tin_Body *island = tin_island_find(body);
			if (body != island) {
				render_push_vertex(body->transform.translation);
				render_push_vertex(island->transform.translation);
			}
		}
		render_draw_lines(TIN_VEC3(0.0, 0.0, 1.0));
#endif

#if 0
		render_start_overlay();
		for (int o1 = 0; o1 < num_objects; o1++) {
			for (int o2 = o1+1; o2 < num_objects; o2++) {
				void *payload;
				if (!tin_find_pair(&scene.arbiters, (uintptr_t)objects[o1].body, (uintptr_t)objects[o2].body, &payload)) {
					continue;
				}
				Tin_Arbiter *arbiter = payload;
				for (int c = 0; c < arbiter->numContacts; c++) {
					Tin_Contact *contact = &arbiter->contacts[c];
					if (contact->separation > 0.0) continue;
					render_push_vertex(tin_fwtrf_point(&arbiter->body1->transform, contact->rel1));
					render_push_vertex(tin_fwtrf_point(&arbiter->body2->transform, contact->rel2));
				}
			}
		}
		render_draw_lines(TIN_VEC3(1.0, 0.0, 0.0));

		render_start_overlay();
		for (int o1 = 0; o1 < num_objects; o1++) {
			for (int o2 = o1+1; o2 < num_objects; o2++) {
				void *payload;
				if (!tin_find_pair(&scene.arbiters, (uintptr_t)objects[o1].body, (uintptr_t)objects[o2].body, &payload)) {
					continue;
				}
				Tin_Arbiter *arbiter = payload;
				for (int c = 0; c < arbiter->numContacts; c++) {
					Tin_Contact *contact = &arbiter->contacts[c];
					if (!(contact->separation > 0.0)) continue;
					render_push_vertex(tin_fwtrf_point(&arbiter->body1->transform, contact->rel1));
					render_push_vertex(tin_fwtrf_point(&arbiter->body2->transform, contact->rel2));
				}
			}
		}
		render_draw_lines(TIN_VEC3(0.0, 1.0, 0.0));

		render_start_overlay();
		for (int o1 = 0; o1 < num_objects; o1++) {
			for (int o2 = o1+1; o2 < num_objects; o2++) {
				void *payload;
				if (!tin_find_pair(&scene.arbiters, (uintptr_t)objects[o1].body, (uintptr_t)objects[o2].body, &payload)) {
					continue;
				}
				Tin_Arbiter *arbiter = payload;
				for (int c = 0; c < arbiter->numContacts; c++) {
					Tin_Contact *contact = &arbiter->contacts[c];
					if (contact->separation > 0.0) continue;

					const Tin_Scalar friction = 0.2;

					Tin_Vec3 r1 = tin_fwtrf_dir(&arbiter->body1->transform,
							tin_scale_v3(arbiter->body1->transform.scale, contact->rel1));
					Tin_Vec3 r2 = tin_fwtrf_dir(&arbiter->body2->transform,
							tin_scale_v3(arbiter->body2->transform.scale, contact->rel2));

					//Tin_Vec3 relVel1 = tin_add_v3(arbiter->body1->velocity, tin_cross_v3(arbiter->body1->angularVelocity, r1));
					//Tin_Vec3 relVel2 = tin_add_v3(arbiter->body2->velocity, tin_cross_v3(arbiter->body2->angularVelocity, r2));
					//Tin_Vec3 relVel = tin_sub_v3(relVel1, relVel2);
					Tin_Vec3 relVel = tin_sub_v3(arbiter->body1->velocity, arbiter->body2->velocity);

					Tin_Vec3 frictionDir = tin_gram_schmidt(contact->normal, relVel);
					if (tin_dot_v3(frictionDir, frictionDir) < 0.001) {
						frictionDir = tin_cross_v3(contact->normal, TIN_VEC3(1.0, 0.0, 0.0));
						if (tin_dot_v3(frictionDir, frictionDir) == 0.0) {
							frictionDir = tin_cross_v3(contact->normal, TIN_VEC3(0.0, 1.0, 0.0));
						}
					}
					frictionDir = tin_normalize_v3(frictionDir);

					Tin_Vec3 orthoDir = tin_cross_v3(contact->normal, frictionDir);

					Tin_Vec3 position = tin_scale_v3(0.5, tin_add_v3(tin_add_v3(arbiter->body1->transform.translation, r1), tin_add_v3(arbiter->body2->transform.translation, r2)));
					render_push_vertex(position);
					render_push_vertex(tin_saxpy_v3(10.0 * contact->ineqAccum * friction, frictionDir, position));
					render_push_vertex(position);
					render_push_vertex(tin_saxpy_v3(10.0 * contact->ineqAccum * friction, orthoDir, position));
				}
			}
		}
		render_draw_lines(TIN_VEC3(0.5, 0.5, 0.5));
#endif

		render_proj_matrix = mat4_orthographic(0, width, 0, height, -1.0, 1.0);
		render_view_matrix = mat4_from_transform(&ident);
		render_start_overlay();

		/* Draw Crosshair */
		int cx = width / 2;
		int cy = height / 2;
		render_push_vertex((Tin_Vec3) {{ cx-6, cy }});
		render_push_vertex((Tin_Vec3) {{ cx+5, cy }});
		render_push_vertex((Tin_Vec3) {{ cx, cy-5 }});
		render_push_vertex((Tin_Vec3) {{ cx, cy+6 }});
		render_draw_lines((Tin_Vec3) {{ 1.0f, 1.0f, 1.0f }});

		/* Draw timing bar graph */
		const Tin_Scalar barY1 = height - 20;
		const Tin_Scalar barY2 = height;
		const Tin_Vec3 barColors[] = {
			{{ 1.0, 0.0, 0.0 }},
			{{ 0.5, 0.5, 0.0 }},
			{{ 0.0, 1.0, 0.0 }},
			{{ 0.0, 0.5, 0.5 }},
			{{ 0.0, 0.0, 1.0 }},
			{{ 0.5, 0.0, 0.5 }},
		};
		double timingTotal = 0.0;
		for (int t = 0; t < 6; t++) {
			timingTotal += accumTimings[t];
		}
		if (timingTotal) {
			Tin_Scalar barX = 0.0;
			for (int t = 0; t < 6; t++) {
				Tin_Scalar barLength = accumTimings[t] / timingTotal * width;
				render_push_vertex(TIN_VEC3(barX, barY1, 0));
				render_push_vertex(TIN_VEC3(barX, barY2, 0));
				render_push_vertex(TIN_VEC3(barX+barLength, barY1, 0));
				render_push_vertex(TIN_VEC3(barX+barLength, barY2, 0));
				render_draw_triangle_strip(barColors[t]);
				barX += barLength;
			}
		}

		glfwSwapBuffers(window);
	}

	render_deinit();
	glfwTerminate();
	return 0;
}

