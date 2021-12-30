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
	tin_body body;
	Model   *model;
} Object;

static Camera camera;
static Object objects[MAX_OBJECTS];
static int    num_objects;

static tin_polytope cone_polytope;
static Model        cone_model;

static void
init_cone(int tessel)
{
	int nverts   = tessel + 2;
	int nindices = 6 * tessel;
	tin_vec3 *verts   = calloc(nverts,   sizeof (tin_vec3));
	GLushort *indices = calloc(nindices, sizeof (GLushort));

	for (int v = 0; v < tessel; v++) {
		float angle = 2.0f * M_PI * v / tessel;
		verts[v] = (tin_vec3) {{ cosf(angle), -1.0f, sinf(angle) }};
	}
	verts[tessel + 0] = (tin_vec3) {{ 0.0f, 1.0f, 0.0f }};
	verts[tessel + 1] = (tin_vec3) {{ 0.0f, -1.0f, 0.0f }};

	for (int v = 0; v < tessel; v++) {
		int w = (v + 1) % tessel;
		indices[6*v+0] = v;
		indices[6*v+1] = w;
		indices[6*v+2] = tessel + 0;
		indices[6*v+3] = w;
		indices[6*v+4] = v;
		indices[6*v+5] = tessel + 1;
	}

	cone_polytope.vertices     = verts;
	cone_polytope.num_vertices = nverts;

	cone_model = render_make_model(nverts, verts, nindices, indices);
	
	free(indices);
}

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

	case GLFW_KEY_W: camera.f = action == GLFW_PRESS; break;
	case GLFW_KEY_A: camera.l = action == GLFW_PRESS; break;
	case GLFW_KEY_S: camera.b = action == GLFW_PRESS; break;
	case GLFW_KEY_D: camera.r = action == GLFW_PRESS; break;
	
	case GLFW_KEY_LEFT_ALT:
		if (action == GLFW_RELEASE) {
			int mode = glfwGetInputMode(window, GLFW_CURSOR);
			if (mode == GLFW_CURSOR_NORMAL) {
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
				glfwGetCursorPos(window, &camera.cursor_x, &camera.cursor_y);
			} else {
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			}
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

        float dx = x - camera.cursor_x;
        float dy = y - camera.cursor_y;

        camera.yaw -= dx * 0.003f;
        while (camera.yaw < 0.0f) {
                camera.yaw += 2.0f * M_PI;
        }
        while (camera.yaw >= 2.0f * M_PI) {
                camera.yaw -= 2.0f * M_PI;
        }

        camera.pitch -= dy * 0.003f;
        if (camera.pitch < -0.5f * M_PI) {
                camera.pitch = -0.5f * M_PI;
        }
        if (camera.pitch > 0.5f * M_PI) {
                camera.pitch = 0.5f * M_PI;
        }

        camera.cursor_x = x;
        camera.cursor_y = y;
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
	gladLoadGL(glfwGetProcAddress);
	glEnable(GL_DEPTH_TEST);
	render_init();
	init_cone(16);

	camera.quat = (tin_quat) { {{ 0.0f, 0.0f, 0.0f }}, 1.0f };
	glfwGetCursorPos(window, &camera.cursor_x, &camera.cursor_y);
	camera.position.c[2] = 5.0f;

	objects[num_objects++] = (Object) {
		{
			{
				tin_make_qt((tin_vec3) {{ 0.0f, 1.0f, 0.0f }}, 0.0f),
				(tin_vec3) {{ 0.0f, 0.0f, 0.0f }},
				1.0f
			},
			&cone_polytope,
			TIN_CONVEX
		},
		&cone_model
	};

	objects[num_objects++] = (Object) {
		{
			{
				tin_make_qt((tin_vec3) {{ 0.0f, 1.0f, 0.0f }}, 0.0f),
				(tin_vec3) {{ 0.75f, 1.0f, 0.0f }},
				0.75f
			},
			&cone_polytope,
			TIN_CONVEX
		},
		&cone_model
	};

	tin_portal portal;
	{
		tin_polysum sum = {
			objects[0].body.shape_params,
			&objects[0].body.transform,
			objects[1].body.shape_params,
			&objects[1].body.transform,
		};
		tin_ray ray;
		ray.origin = tin_sub_v3(
				objects[0].body.transform.translation,
				objects[1].body.transform.translation);
		ray.dir = tin_normalize_v3(tin_neg_v3(ray.origin));

		if (tin_construct_portal(&sum, &ray, &portal)) {
			printf("Got portal\n");
			tin_refine_portal(&sum, &ray, &portal);
		} else {
			printf("Got no portal\n");
		}
	}

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();
		camera_update(&camera, 1.0f / 60.0f);

		int width, height;
		glfwGetWindowSize(window, &width, &height);
		glViewport(0, 0, width, height);

		glClearColor(0.7f, 0.9f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		tin_transform camtrf;
		camera_transform(&camera, &camtrf);

		tin_transform ident = { .rotation = { {{ 0.0f, 0.0f, 0.0f }}, 1.0f }, .scale = 1.0f };

		render_start_models();
		render_proj_matrix = mat4_perspective(70.0f, (float) width / height, 1.0f, 100.0f);
		render_view_matrix = mat4_from_inverse_transform(&camtrf);

		for (int o = 0; o < num_objects; o++) {
			const Object *obj = &objects[o];
			tin_polytope dot;
			dot.num_vertices = 1;
			dot.vertices = malloc(1 * sizeof *dot.vertices);
			dot.vertices[0] = (tin_vec3) {{ 0.0f, 0.0f, 0.0f }};
			tin_polysum sum;
			sum.former_polytope = obj->body.shape_params;
			sum.former_transform = &obj->body.transform;
			sum.latter_polytope = &dot;
			sum.latter_transform = &ident;
			tin_ray ray = { .origin = camtrf.translation, .dir = tin_apply_qt(camtrf.rotation, (tin_vec3) {{ 0.0f, 0.0f, -1.0f }}) };
			tin_vec3 color;
			tin_portal temp_portal;
			if (tin_construct_portal(&sum, &ray, &temp_portal)) {
				color = (tin_vec3) {{ 1.0f, 0.0f, 0.0f }};
			} else {
				color = (tin_vec3) {{ 1.0f, 1.0f, 1.0f }};
			}
			render_draw_model(obj->model, &obj->body.transform, color);
		}

		tin_vec3 color = {{ 0.0f, 1.0f, 1.0f }};
		tin_transform trf = ident;
		trf.scale = 0.1f;
		trf.translation = tin_fwtrf_point(&objects[0].body.transform, portal.a.rel_former);
		render_draw_model(&cone_model, &trf, color);
		trf.translation = tin_fwtrf_point(&objects[0].body.transform, portal.b.rel_former);
		render_draw_model(&cone_model, &trf, color);
		trf.translation = tin_fwtrf_point(&objects[0].body.transform, portal.c.rel_former);
		render_draw_model(&cone_model, &trf, color);
		trf.translation = tin_fwtrf_point(&objects[1].body.transform, portal.a.rel_latter);
		render_draw_model(&cone_model, &trf, color);
		trf.translation = tin_fwtrf_point(&objects[1].body.transform, portal.b.rel_latter);
		render_draw_model(&cone_model, &trf, color);
		trf.translation = tin_fwtrf_point(&objects[1].body.transform, portal.c.rel_latter);
		render_draw_model(&cone_model, &trf, color);

		render_start_overlay();

		render_proj_matrix = mat4_orthographic(width, height);
		render_view_matrix = mat4_from_transform(&ident);

		int cx = width / 2;
		int cy = height / 2;
		render_push_vertex((tin_vec3) {{ cx-6, cy }});
		render_push_vertex((tin_vec3) {{ cx+5, cy }});
		render_push_vertex((tin_vec3) {{ cx, cy-5 }});
		render_push_vertex((tin_vec3) {{ cx, cy+6 }});
		render_draw_lines((tin_vec3) {{ 1.0f, 1.0f, 1.0f }});

		glfwSwapBuffers(window);
	}

	render_deinit();
	glfwTerminate();
	return 0;
}

