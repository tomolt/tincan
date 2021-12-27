typedef struct {
	GLsizei    index_count;
	GLsizeiptr base_index;
	GLint      base_vertex;
} Model;

void  render_init  (void);
void  render_deinit(void);
Model render_make_model(int nverts, const tin_vec3 *verts, int nindices, const GLushort *indices);
void  render_model(const Model *model, const tin_transform *transform, const tin_transform *camera_transform, int width, int height);

