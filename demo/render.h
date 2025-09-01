typedef struct {
	GLsizei    index_count;
	GLsizeiptr base_index;
	GLint      base_vertex;
} Model;

extern Mat4 render_view_matrix;
extern Mat4 render_proj_matrix;

void render_init(void);
void render_deinit(void);

Model render_make_model(int nverts, const Tin_Vec3 *verts, int nindices, const GLushort *indices);

void render_start_models(void);
void render_draw_model(const Model *model, const Tin_Transform *transform, Tin_Vec3 color);

void render_start_overlay(void);
void render_push_vertex(Tin_Vec3 v);
void render_draw_lines(Tin_Vec3 color);
void render_draw_triangle_strip(Tin_Vec3 color);

