#include <stdlib.h>
#include <stdio.h>

#include <glad/gl.h>

#include <tincan.h>

#include "render.h"

static GLuint shader_prog;

static const char *shader_vert_src =
	"#version 330 core\n"
	"uniform mat4 mvp_matrix;\n"
	"uniform mat4 normal_matrix;\n"
	"layout(location=0) in vec3 position;\n"
	"layout(location=1) in vec3 normal;\n"
	"out vec3 frag_normal;\n"
	"void main() {\n"
	"	gl_Position = mvp_matrix * vec4(position, 1.0);\n"
	"	frag_normal = vec3(normal_matrix * vec4(normal, 0.0));\n"
	"}\n";

static const char *shader_frag_src =
	"#version 330 core\n"
	"uniform vec4 base_color;\n"
	"uniform vec3 light_dir;\n"
	"in vec3 frag_normal;\n"
	"out vec4 frag_color;\n"
	"void main() {\n"
	"	float diffuse = max(dot(frag_normal, -light_dir), 0.0);\n"
	"	frag_color = base_color * (0.8 * diffuse + 0.2);\n"
	"}\n";

GLuint
shader_compile(const char *source, GLenum type)
{
	GLuint shader = glCreateShader(type);
	glShaderSource(shader, 1, (const GLchar * const *) &source, NULL);
	glCompileShader(shader);

	GLint success;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (!success) {
		char msg_buf[1024];
		GLsizei msg_len;
		glGetShaderInfoLog(shader, sizeof msg_buf, &msg_len, msg_buf);
		fprintf(stderr, "glsl compile error: %s\n", msg_buf);
	}

	return shader;
}

GLuint
shader_link(GLuint vert, GLuint frag)
{
	GLuint prog = glCreateProgram();
	
	glAttachShader(prog, vert);
	glAttachShader(prog, frag);
	
	glLinkProgram(prog);

	GLint success;
	glGetProgramiv(prog, GL_LINK_STATUS, &success);
	if (!success) {
		char msg_buf[1024];
		GLsizei msg_len;
		glGetProgramInfoLog(prog, sizeof msg_buf, &msg_len, msg_buf);
		fprintf(stderr, "glsl link error: %s\n", msg_buf);
	}
	
	glDetachShader(prog, vert);
	glDetachShader(prog, frag);
	
	return prog;
}

void
render_init(void)
{
	GLuint vert = shader_compile(shader_vert_src, GL_VERTEX_SHADER);
	GLuint frag = shader_compile(shader_frag_src, GL_FRAGMENT_SHADER);

	shader_prog = shader_link(vert, frag);

	glDeleteShader(vert);
	glDeleteShader(frag);
}

void
render_deinit(void)
{
	glDeleteProgram(shader_prog);
}

void
render_model(ModelRef model, const tin_transform *transform)
{
}

