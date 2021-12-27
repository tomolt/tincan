#include <tincan.h>

#include <stdlib.h>
#include <stdio.h>

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include "render.h"

const GLuint WIDTH = 800, HEIGHT = 600;

void
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

	GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, "tincan physics demo", NULL, NULL);
	glfwSetKeyCallback(window, key_callback);

	glfwMakeContextCurrent(window);
	gladLoadGL(glfwGetProcAddress);
	render_init();

	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		glClearColor(0.7f, 0.9f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);



		glfwSwapBuffers(window);
	}

	render_deinit();
	glfwTerminate();
	return 0;
}

