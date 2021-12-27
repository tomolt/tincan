#include <stdlib.h>
#include <stdio.h>

#include <glad/gl.h>

#include <tincan.h>

#include "mat4.h"
#include "render.h"

#define VBO_CAPACITY 4096
#define IBO_CAPACITY 4096

static GLuint model_prog;
static GLuint model_vbo;
static GLuint model_ibo;
static GLuint model_vao;
static GLuint model_uniforms[4];

static GLsizei model_vbo_count;
static GLsizei model_ibo_count;

static const char *model_vert_src =
	"#version 330 core\n"
	"uniform mat4 model_view_matrix;\n"
	"uniform mat4 proj_matrix;\n"
	"layout(location=0) in vec3 position;\n"
	"out vec4 view_position;\n"
	"void main() {\n"
	"	view_position = model_view_matrix * vec4(position, 1.0);\n"
	"	gl_Position = proj_matrix * view_position;\n"
	"}\n";

static const char *model_frag_src =
	"#version 330 core\n"
	"uniform mat4 model_view_matrix;\n"
	"uniform vec3 light_dir;\n"
	"uniform vec4 base_color;\n"
	"in vec4 view_position;\n"
	"out vec4 frag_color;\n"
	"void main() {\n"
	"	vec3 tangent = dFdx(view_position.xyz);\n"
	"	vec3 bitangent = dFdy(view_position.xyz);\n"
	"	vec3 normal = normalize(cross(tangent, bitangent));\n"
	"	vec3 view_light_dir = vec3(model_view_matrix * vec4(light_dir, 0.0));\n"
	"	float diffuse = max(dot(normal, -view_light_dir), 0.0);\n"
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
	GLuint vert = shader_compile(model_vert_src, GL_VERTEX_SHADER);
	GLuint frag = shader_compile(model_frag_src, GL_FRAGMENT_SHADER);
	model_prog = shader_link(vert, frag);
	glDeleteShader(vert);
	glDeleteShader(frag);
	model_uniforms[0] = glGetUniformLocation(model_prog, "model_view_matrix");
	model_uniforms[1] = glGetUniformLocation(model_prog, "proj_matrix");
	model_uniforms[2] = glGetUniformLocation(model_prog, "light_dir");
	model_uniforms[3] = glGetUniformLocation(model_prog, "base_color");

	glGenBuffers(1, &model_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
	glBufferData(GL_ARRAY_BUFFER, VBO_CAPACITY, NULL, GL_STATIC_DRAW);
	model_vbo_count = 0;

	glGenBuffers(1, &model_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IBO_CAPACITY, NULL, GL_STATIC_DRAW);
	model_ibo_count = 0;

	glGenVertexArrays(1, &model_vao);
	glBindVertexArray(model_vao);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_TRUE, 0, (void *) 0);
}

void
render_deinit(void)
{
	glDeleteVertexArrays(1, &model_vao);
	glDeleteBuffers(1, &model_vbo);
	glDeleteBuffers(1, &model_ibo);
	glDeleteProgram(model_prog);
}

Model
render_make_model(int nverts, const tin_vec3 *verts, int nindices, const GLushort *indices)
{
	Model model = {
		.base_index  = model_ibo_count,
		.base_vertex = model_vbo_count,
		.index_count = nindices
	};

	glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
	glBufferSubData(GL_ARRAY_BUFFER,
		model_vbo_count * sizeof (tin_vec3),
		nverts * sizeof (tin_vec3), verts);
	model_vbo_count += nverts;

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model_ibo);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,
		model_ibo_count * sizeof (GLushort),
		nindices * sizeof (GLushort), indices);
	model_ibo_count += nindices;

	return model;
}

void
render_model(const Model *model, const tin_transform *transform, const tin_transform *camera_transform, int width, int height)
{
	glUseProgram(model_prog);
	glBindVertexArray(model_vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model_ibo);

	Mat4 model_matrix = mat4_from_transform(transform);
	Mat4 view_matrix = mat4_from_inverse_transform(camera_transform);
	Mat4 model_view_matrix = mat4_multiply(view_matrix, model_matrix);
	glUniformMatrix4fv(model_uniforms[0], 1, GL_FALSE, model_view_matrix.c);
	Mat4 proj_matrix = mat4_projection(70.0f, (float) width / height, 1.0f, 100.0f);
	glUniformMatrix4fv(model_uniforms[1], 1, GL_FALSE, proj_matrix.c);
	tin_vec3 light_dir = tin_normalize_v3((tin_vec3) {{ -1.0f, -1.0f, 0.0f }});
	glUniform3fv(model_uniforms[2], 1, light_dir.c);
	glUniform4f(model_uniforms[3], 1.0f, 1.0f, 1.0f, 1.0f);

	glDrawElementsBaseVertex(GL_TRIANGLES, model->index_count,
		GL_UNSIGNED_SHORT, (void *) (model->base_index * sizeof (GLushort)), model->base_vertex);
}

