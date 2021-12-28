#include <stdlib.h>
#include <stdio.h>

#include <glad/gl.h>

#include <tincan.h>

#include "mat4.h"
#include "render.h"

#define MODEL_VBO_CAPACITY   4096
#define MODEL_IBO_CAPACITY   4096
#define OVERLAY_VBO_CAPACITY 4096

static GLuint  model_prog;
static GLuint  model_uniforms[4];
static GLuint  model_vbo;
static GLuint  model_ibo;
static GLuint  model_vao;
static GLsizei model_vbo_count;
static GLsizei model_ibo_count;

static GLuint  overlay_prog;
static GLuint  overlay_uniforms[3];
static GLuint  overlay_vbo;
static GLuint  overlay_vao;
static GLsizei overlay_vbo_start;
static GLsizei overlay_vbo_count;

Mat4 render_view_matrix;
Mat4 render_proj_matrix;

static const char *common_vert_src =
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

static const char *overlay_frag_src =
	"#version 330 core\n"
	"uniform vec4 base_color;\n"
	"in vec4 view_position;\n"
	"out vec4 frag_color;\n"
	"void main() {\n"
	"	frag_color = base_color;\n"
	"}\n";

static GLuint
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

static GLuint
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

static void
init_progs(void)
{
	GLuint common_vert  = shader_compile(common_vert_src,  GL_VERTEX_SHADER);
	GLuint model_frag   = shader_compile(model_frag_src,   GL_FRAGMENT_SHADER);
	GLuint overlay_frag = shader_compile(overlay_frag_src, GL_FRAGMENT_SHADER);

	model_prog = shader_link(common_vert, model_frag);
	model_uniforms[0] = glGetUniformLocation(model_prog, "model_view_matrix");
	model_uniforms[1] = glGetUniformLocation(model_prog, "proj_matrix");
	model_uniforms[2] = glGetUniformLocation(model_prog, "light_dir");
	model_uniforms[3] = glGetUniformLocation(model_prog, "base_color");

	overlay_prog = shader_link(common_vert, overlay_frag);
	overlay_uniforms[0] = glGetUniformLocation(overlay_prog, "model_view_matrix");
	overlay_uniforms[1] = glGetUniformLocation(overlay_prog, "proj_matrix");
	overlay_uniforms[2] = glGetUniformLocation(overlay_prog, "base_color");
	
	glDeleteShader(common_vert);
	glDeleteShader(model_frag);
	glDeleteShader(overlay_frag);
}

void
render_init(void)
{
	init_progs();

	/* model */

	glGenBuffers(1, &model_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
	glBufferData(GL_ARRAY_BUFFER, MODEL_VBO_CAPACITY, NULL, GL_STATIC_DRAW);
	model_vbo_count = 0;

	glGenBuffers(1, &model_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, MODEL_IBO_CAPACITY, NULL, GL_STATIC_DRAW);
	model_ibo_count = 0;

	glGenVertexArrays(1, &model_vao);
	glBindVertexArray(model_vao);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_TRUE, 0, (void *) 0);

	/* overlay */

	glGenBuffers(1, &overlay_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, overlay_vbo);

	glGenVertexArrays(1, &overlay_vao);
	glBindVertexArray(overlay_vao);
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

	glDeleteVertexArrays(1, &overlay_vao);
	glDeleteBuffers(1, &overlay_vbo);
	glDeleteProgram(overlay_prog);
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
render_start_models(void)
{
	glUseProgram(model_prog);
	glBindVertexArray(model_vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model_ibo);
	glEnable(GL_DEPTH_TEST);

	tin_vec3 light_dir = tin_normalize_v3((tin_vec3) {{ -1.0f, -1.0f, 0.0f }});
	glUniform3fv(model_uniforms[2], 1, light_dir.c);
}

void
render_draw_model(const Model *model, const tin_transform *transform, tin_vec3 color)
{
	Mat4 model_matrix = mat4_from_transform(transform);
	Mat4 model_view_matrix = mat4_multiply(render_view_matrix, model_matrix);
	glUniformMatrix4fv(model_uniforms[0], 1, GL_FALSE, model_view_matrix.c);
	glUniformMatrix4fv(model_uniforms[1], 1, GL_FALSE, render_proj_matrix.c);
	glUniform4f(model_uniforms[3], color.c[0], color.c[1], color.c[2], 1.0f);

	glDrawElementsBaseVertex(GL_TRIANGLES, model->index_count,
		GL_UNSIGNED_SHORT, (void *) (model->base_index * sizeof (GLushort)), model->base_vertex);
}

void
render_start_overlay(void)
{
	glUseProgram(overlay_prog);
	glBindVertexArray(overlay_vao);
	glBindBuffer(GL_ARRAY_BUFFER, overlay_vbo);
	glBufferData(GL_ARRAY_BUFFER, OVERLAY_VBO_CAPACITY, NULL, GL_STREAM_DRAW);
	glDisable(GL_DEPTH_TEST);
	overlay_vbo_start = 0;
	overlay_vbo_count = 0;
}

void
render_push_vertex(tin_vec3 v)
{
	glBufferSubData(GL_ARRAY_BUFFER,
		(overlay_vbo_start + overlay_vbo_count) * sizeof (tin_vec3),
		sizeof (tin_vec3), v.c);
	overlay_vbo_count++;
}

void
render_draw_lines(tin_vec3 color)
{
	glUniformMatrix4fv(overlay_uniforms[0], 1, GL_FALSE, render_view_matrix.c);
	glUniformMatrix4fv(overlay_uniforms[1], 1, GL_FALSE, render_proj_matrix.c);
	glUniform4f(overlay_uniforms[2], color.c[0], color.c[1], color.c[2], 1.0f);
	glDrawArrays(GL_LINES, overlay_vbo_start, overlay_vbo_count);
	overlay_vbo_start += overlay_vbo_count;
	overlay_vbo_count  = 0;
}

