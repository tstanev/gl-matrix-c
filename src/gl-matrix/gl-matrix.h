#pragma once
#include "defs.h"

/*
#include "vec2.c"
#include "vec3.c"
*/


//#include "vec4.c"
vec4_p vec4_normalize(vec4_p out, vec4_cp a);
vec4_p vec4_scale(vec4_p out, vec4_cp a, scalar_t b);

/*
#include "mat2.c"
#include "mat2d.c"
#include "mat3.c"
#include "mat4.c"
#include "quat.c"
#include "quat2.c"
*/

//#include "utils.c"
scalar_t random1();
vec4_p vec4_random(vec4_p out);