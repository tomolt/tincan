typedef struct {
	GLsizei  index_count;
	intptr_t base_index;
	GLint    base_vertex;
} Model;

void render_init  (void);
void render_deinit(void);
void render_model(const Model *model, const tin_transform *transform, int width, int height);
