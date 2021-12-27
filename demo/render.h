typedef struct {
} ModelRef;

void render_init  (void);
void render_deinit(void);
void render_model (ModelRef model, const tin_transform *transform);
