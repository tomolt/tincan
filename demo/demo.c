#include <tincan.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
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

static Camera camera;
static Object objects[MAX_OBJECTS];
static int    num_objects;

static Tin_Shape cone_shape;
static Model     cone_model;
static Tin_Shape cube_shape;
static Model     cube_model;

static float time_multiplier = 1.0f;

static void
init_cone(int tessel)
{
	int nverts   = tessel + 2;
	int nindices = 6 * tessel;
	Tin_Vec3 *verts   = calloc(nverts,   sizeof (Tin_Vec3));
	GLushort *indices = calloc(nindices, sizeof (GLushort));

	for (int v = 0; v < tessel; v++) {
		float angle = 2.0f * M_PI * v / tessel;
		verts[v] = (Tin_Vec3) {{ cosf(angle), 0.5f, sinf(angle) }};
	}
	verts[tessel + 0] = (Tin_Vec3) {{ 0.0f, -1.5f, 0.0f }};
	verts[tessel + 1] = (Tin_Vec3) {{ 0.0f, 0.5f, 0.0f }};

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
		0, 1, 2, 1, 2, 3,
		4, 5, 6, 5, 6, 7,
		0, 1, 4, 1, 4, 5,
		2, 3, 6, 3, 6, 7,
		0, 2, 4, 2, 4, 6,
		1, 3, 5, 3, 5, 7,
	};

	cube_shape.kind = TIN_POLYTOPE;
	cube_shape.polytope.vertices = verts;
	cube_shape.polytope.numVertices = 8;
	cube_shape.radius = sqrtf(3.0f);
	cube_shape.invInertia = (Tin_Vec3) {{ 6.0f, 6.0f, 6.0f }};
	
	cube_model = render_make_model(8, verts, 6 * 6, indices);
}

static Camera *cam = &camera;

static void
key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	(void) scancode, (void) mods;

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
			Tin_Vec3 forward = tin_apply_qt(cam->quat, (Tin_Vec3) {{ 0.0f, 0.0f, -1.0f }});
			Tin_Body *body = tin_add_body(&scene);
			*body = (Tin_Body) {
				body->node,
				{
					cam->quat,
					tin_saxpy_v3(1.0f, forward, cam->position),
					0.2f,
				},
				&cube_shape,
				1.0f / 1.0f,
				tin_scale_v3(5.0f, forward),
				{{ 0.0f, 0.0f, 0.0f }},
			};
			objects[num_objects++] = (Object) {
				body,
				&cube_model,
				{{ 0.5f, 0.5f, 1.0f }}
			};
		}
		break;
	}
}

static void
mouse_callback(GLFWwindow *window, double x, double y)
{
	(void) window;

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

void *
custom_alloc(void *userPointer)
{
	return malloc((size_t) userPointer);
}

void
custom_free(void *userPointer, void *memoryPointer)
{
	(void) userPointer;
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
	init_cone(16);
	init_cube();

	camera.quat = (Tin_Quat) { {{ 0.0f, 0.0f, 0.0f }}, 1.0f };
	glfwGetCursorPos(window, &camera.cursor_x, &camera.cursor_y);
	camera.position.c[2] = 5.0f;

	TIN_LIST_INIT(scene.bodies);
	TIN_LIST_INIT(scene.arbiters);
	TIN_LIST_INIT(scene.joints);
	scene.bodyAllocator = (Tin_Allocator) { (void *) sizeof (Tin_Body), custom_alloc, custom_free };
	scene.arbiterAllocator = (Tin_Allocator) { (void *) sizeof (Tin_Arbiter), custom_alloc, custom_free };
	scene.jointAllocator = (Tin_Allocator) { (void *) sizeof (Tin_Joint), custom_alloc, custom_free };

	Tin_Body *body1 = tin_add_body(&scene);
	*body1 = (Tin_Body) {
		body1->node,
		{
			tin_make_qt((Tin_Vec3) {{ 0.0f, 1.0f, 0.0f }}, 0.0f),
			(Tin_Vec3) {{ 0.0f, 1.5f, 0.0f }},
			1.0f
		},
		&cone_shape,
		1.0f / 1.0f,
		//{{ 0.0f, 0.0f, -0.5f }},
		{{ 0.0f, 0.0f, 0.0f }},
		{{ 0.0f, 0.0f, 0.0f }},
	};
	objects[num_objects++] = (Object) {
		body1,
		&cone_model,
		{{ 0.5f, 1.0f, 0.5f }}
	};

	Tin_Body *body2 = tin_add_body(&scene);
	*body2 = (Tin_Body) {
		body2->node,
		{
			tin_make_qt((Tin_Vec3) {{ 0.0f, 1.0f, 0.0f }}, 0.0f),
			(Tin_Vec3) {{ 0.0f, -1.0f, 0.0f }},
			1.0f
		},
		&cube_shape,
		1.0f / 3.0f,
		{{ -3.0f, 1.0f, 0.0f }},
		{{ 0.0f, 0.0f, 0.0f }},
		//{{ 0.1f, 0.1f, 0.0f }},
	};
	objects[num_objects++] = (Object) {
		body2,
		&cube_model,
		{{ 1.0f, 0.5f, 0.5f }}
	};

	Tin_Body *body3 = tin_add_body(&scene);
	*body3 = (Tin_Body) {
		body3->node,
		{
			tin_make_qt((Tin_Vec3) {{ 0.0f, 1.0f, 0.0f }}, 0.0f),
			(Tin_Vec3) {{ 0.0f, -23.0f, 0.0f }},
			20.0f
		},
		&cube_shape,
		0.0f,
		{{ 0.0f, 0.0f, 0.0f }},
		{{ 0.0f, 0.0f, 0.0f }},
	};
	objects[num_objects++] = (Object) {
		body3,
		&cube_model,
		{{ 1.0f, 1.0f, 1.0f }}
	};

	Tin_Joint *joint = tin_add_joint(&scene);
	joint->body1 = body1;
	joint->body2 = body2;
	joint->relTo1 = (Tin_Vec3){{ 0.0f, -1.505f, 0.0f }};
	joint->relTo2 = (Tin_Vec3){{ 0.0f, 1.005f, 0.0f }};

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		Tin_Scalar dt = 1.0f / 60.0f * time_multiplier;

		camera_update(&camera, 1.0f / 60.0f);

		int width, height;
		glfwGetWindowSize(window, &width, &height);
		glViewport(0, 0, width, height);

		glClearColor(0.7f, 0.9f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		Tin_Transform camtrf;
		camera_transform(&camera, &camtrf);

		Tin_Transform ident = { .rotation = { {{ 0.0f, 0.0f, 0.0f }}, 1.0f }, .scale = 1.0f };

		render_start_models();
		render_proj_matrix = mat4_perspective(70.0f, (float) width / height, 1.0f, 100.0f);
		render_view_matrix = mat4_from_inverse_transform(&camtrf);

		tin_simulate(&scene, dt);

		/* draw contact points */
#if 0
		for (int a = 0; a < num_objects; a++) {
			for (int b = a + 1; b < num_objects; b++) {
				Tin_Arbiter *arbiter = &arbiters[a * MAX_OBJECTS + b];
				for (int i = 0; i < arbiter->numContacts; i++) {
					const Tin_Contact *c = &arbiter->contacts[i];
					Tin_Vec3 color = {{ 0.0f, 1.0f, 1.0f }};
					Tin_Transform trf = ident;
					trf.scale = 0.1f;
					trf.translation = tin_fwtrf_point(&objects[a].body.transform, c->rel1);
					render_draw_model(&cone_model, &trf, color);
					trf.translation = tin_fwtrf_point(&objects[b].body.transform, c->rel2);
					render_draw_model(&cube_model, &trf, color);
				}
			}
		}
#endif

#if 0
		/* draw joints */
		TIN_FOR_EACH(joint, scene.joints, Tin_Joint, node) {
			const Tin_Vec3 color1 = {{ 1.0f, 1.0f, 1.0f }};
			const Tin_Vec3 color2 = {{ 0.0f, 0.0f, 0.0f }};

			Tin_Transform trf = ident;
			trf.scale = 0.1f;
			
			trf.translation = tin_fwtrf_point(&joint->body1->transform, joint->relTo1);
			render_draw_model(&cube_model, &trf, color1);
			trf.translation = tin_sub_v3(trf.translation, joint->debugImpulse);
			render_draw_model(&cone_model, &trf, color1);

			trf.translation = tin_fwtrf_point(&joint->body2->transform, joint->relTo2);
			render_draw_model(&cube_model, &trf, color2);
			trf.translation = tin_add_v3(trf.translation, joint->debugImpulse);
			render_draw_model(&cone_model, &trf, color2);
		}
#endif

		for (int o = 0; o < num_objects; o++) {
			const Object *obj = &objects[o];
			Tin_Polytope dot;
			dot.numVertices = 1;
			dot.vertices = malloc(1 * sizeof *dot.vertices);
			dot.vertices[0] = (Tin_Vec3) {{ 0.0f, 0.0f, 0.0f }};
			/* ray picking (outdated code) */
#if 0
			Tin_Polysum sum;
			sum.polytope1 = obj->body->shape->polytope;
			sum.transform1 = &obj->body.transform;
			sum.polytope2 = &dot;
			sum.transform2 = &ident;
			Tin_Ray ray = { .origin = camtrf.translation, .dir = tin_apply_qt(camtrf.rotation, (Tin_Vec3) {{ 0.0f, 0.0f, -1.0f }}) };
			Tin_Portal temp_portal;
			if (tin_construct_portal(&sum, &ray, &temp_portal)) {
				color = (Tin_Vec3) {{ 1.0f, 0.0f, 0.0f }};
			} else {
				color = (Tin_Vec3) {{ 1.0f, 1.0f, 1.0f }};
			}
#endif
			render_draw_model(obj->model, &obj->body->transform, obj->color);
		}

		render_start_overlay();

		render_proj_matrix = mat4_orthographic(width, height);
		render_view_matrix = mat4_from_transform(&ident);

		int cx = width / 2;
		int cy = height / 2;
		render_push_vertex((Tin_Vec3) {{ cx-6, cy }});
		render_push_vertex((Tin_Vec3) {{ cx+5, cy }});
		render_push_vertex((Tin_Vec3) {{ cx, cy-5 }});
		render_push_vertex((Tin_Vec3) {{ cx, cy+6 }});
		render_draw_lines((Tin_Vec3) {{ 1.0f, 1.0f, 1.0f }});

		glfwSwapBuffers(window);
	}

	render_deinit();
	glfwTerminate();
	return 0;
}

