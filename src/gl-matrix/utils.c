
#include <stdlib.h>
#include "utils.h"

/* Returns a random number in the interval [0, 1) -- very simple implementation, FOR TESTING ONLY */
scalar_t random1() {
	return (scalar_t)((double)rand() / ((double)RAND_MAX+1));
}

/**
 * Generates a random vector with the given scale
 *
 * @param {vec4} out the receiving vector
 * @param {Number} [scale] Length of the resulting vector. If omitted, a unit vector will be returned
 * @returns {vec4} out
 */

vec4_p vec4_random(vec4_p out) {

	//TODO: move this to function params
	scalar_t vectorScale = 1.0;

	out[0] = random1();
	out[1] = random1();
	out[2] = random1();
	out[3] = random1();
	vec4_normalize(out, out);
	vec4_scale(out, out, vectorScale);
	return out;
}
