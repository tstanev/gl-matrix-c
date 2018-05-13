#include "defs.h"

/* "Nice" declarations for convenience of declaring/creating the structures on the stack
 -- not used directly by the gl-matrix API */
typedef struct vec2_t {
	union {
		struct {
			scalar_t x;
			scalar_t y;
		};
		scalar_t p[2];
	};
} vec2_t;

typedef struct vec3_t {
	union {
		struct {
			scalar_t x;
			scalar_t y;
			scalar_t z;
		};
		scalar_t p[3];
	};
} vec3_t;

typedef struct vec4_t {
	union {
		struct {
			scalar_t x;
			scalar_t y;
			scalar_t z;
			scalar_t w;
		};
		scalar_t p[4];
	};
} vec4_t;

typedef struct mat3_t {
	union {
		vec3_t cols[3];
		scalar_t p[9];
	};
} mat3_t;

typedef struct mat4_t {
	union {
		vec4_t cols[4];
		scalar_t p[16];
	};
} mat4_t;

typedef vec4_t quat_t;
