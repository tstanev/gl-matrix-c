#include <stdio.h>
#include <stdlib.h>
#include "gl-matrix/gl-matrix.h"
#include "gl-matrix/nice-types.h"
#include "gl-matrix/utils.h"

int main(int argc, char** argv) {

	srand(17);

	vec3_t v;
	v.x = 1;
	v.y = 1;
	v.z = 1;

	vec4_t r1;
	vec4_t r2;
	vec4_random(r1.p);
	vec4_random(r2.p);

	//double d = vec3_dot(v.p, v.p);
	printf("%d\n", RAND_MAX);
	printf("%lf %lf %lf %lf\n", r1.x, r1.y, r1.z, r1.w);
}