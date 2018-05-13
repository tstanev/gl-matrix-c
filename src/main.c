#include <stdio.h>
#include "gl-matrix/gl-matrix.h"
#include "gl-matrix/nice-types.h"

int main(int argc, char** argv) {

	vec3_t v;
	v.x = 1;
	v.y = 1;
	v.z = 1;

	double d = vec3_dot(v.p, v.p);

	printf("Hello world. %lf\n", d);
}