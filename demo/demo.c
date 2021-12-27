#include <tincan.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include "render.h"

#define MAX_OBJECTS 100

typedef struct {
	tin_body body;
	Model   *model;
} Object;

static Object objects[MAX_OBJECTS];
static int num_objects;

static Model cone_model;

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

	cone_model = render_make_model(nverts, verts, nindices, indices);
	
	free(verts);
	free(indices);
}

static void
key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	(void) scancode, (void) mode;
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
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

	glfwMakeContextCurrent(window);
	gladLoadGL(glfwGetProcAddress);
	glEnable(GL_DEPTH_TEST);
	render_init();
	init_cone(16);

	objects[num_objects++] = (Object) {
		{
			{
				tin_make_qt((tin_vec3) {{ 0.0f, 1.0f, 0.0f }}, 0.0f),
				(tin_vec3) {{ 0.0f, 0.0f, -5.0f }},
				1.0f
			},
			NULL
		},
		&cone_model
	};

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		int width, height;
		glfwGetWindowSize(window, &width, &height);
		glViewport(0, 0, width, height);

		glClearColor(0.7f, 0.9f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		for (int o = 0; o < num_objects; o++) {
			const Object *obj = &objects[o];
			render_model(obj->model, &obj->body.transform, width, height);
		}

		glfwSwapBuffers(window);
	}

	render_deinit();
	glfwTerminate();
	return 0;
}

