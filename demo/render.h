typedef struct {
	GLsizei    index_count;
	GLsizeiptr base_index;
	GLint      base_vertex;
} Model;

extern Mat4 render_view_matrix;
extern Mat4 render_proj_matrix;

void render_init(void);
void render_deinit(void);

Model render_make_model(int nverts, const tin_vec3 *verts, int nindices, const GLushort *indices);

void render_start_models(void);
void render_draw_model(const Model *model, const tin_transform *transform, tin_vec3 color);

void render_start_overlay(void);
void render_push_vertex(tin_vec3 v);
void render_draw_lines(tin_vec3 color);

