/* Copyright (c) 2015, Brandon Jones, Colin MacKenzie IV.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE. */

#include "defs.h"

/**
 * 2 Dimensional Vector
 * @module vec2
 */

/**
 * Copy the values from one vec2 to another
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the source vector
 * @returns {vec2} out
 */
vec2_p vec2_copy(vec2_p out, vec2_cp a) {
	out[0] = a[0];
	out[1] = a[1];
	return out;
}

/**
 * Set the components of a vec2 to the given values
 *
 * @param {vec2} out the receiving vector
 * @param {Number} x X component
 * @param {Number} y Y component
 * @returns {vec2} out
 */
vec2_p vec2_set(vec2_p out, scalar_t x, scalar_t y) {
	out[0] = x;
	out[1] = y;
	return out;
}

/**
 * Adds two vec2's
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec2} out
 */
vec2_p vec2_add(vec2_p out, vec2_cp a, vec2_cp b) {
	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
	return out;
}

/**
 * Subtracts vector b from vector a
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec2} out
 */
vec2_p vec2_sub(vec2_p out, vec2_cp a, vec2_cp b) {
	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
	return out;
}

/**
 * Multiplies two vec2's
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec2} out
 */
vec2_p vec2_mul(vec2_p out, vec2_cp a, vec2_cp b) {
	out[0] = a[0] * b[0];
	out[1] = a[1] * b[1];
	return out;
};

/**
 * Divides two vec2's
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec2} out
 */
vec2_p vec2_div(vec2_p out, vec2_cp a, vec2_cp b) {
	out[0] = a[0] / b[0];
	out[1] = a[1] / b[1];
	return out;
};

/**
 * Math.ceil the components of a vec2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a vector to ceil
 * @returns {vec2} out
 */
vec2_p vec2_ceil(vec2_p out, vec2_cp a) {
	out[0] = ceil(a[0]);
	out[1] = ceil(a[1]);
	return out;
};

/**
 * Math.floor the components of a vec2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a vector to floor
 * @returns {vec2} out
 */
vec2_p vec2_floor(vec2_p out, vec2_cp a) {
	out[0] = floor(a[0]);
	out[1] = floor(a[1]);
	return out;
};

/**
 * Returns the minimum of two vec2's
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec2} out
 */
vec2_p vec2_min(vec2_p out, vec2_cp a, vec2_cp b) {
	out[0] = fmin(a[0], b[0]);
	out[1] = fmin(a[1], b[1]);
	return out;
};

/**
 * Returns the maximum of two vec2's
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec2} out
 */
vec2_p vec2_max(vec2_p out, vec2_cp a, vec2_cp b) {
	out[0] = fmax(a[0], b[0]);
	out[1] = fmax(a[1], b[1]);
	return out;
};

/**
 * Math.round the components of a vec2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a vector to round
 * @returns {vec2} out
 */
vec2_p vec2_round(vec2_p out, vec2_cp a) {
	out[0] = round(a[0]);
	out[1] = round(a[1]);
	return out;
};

/**
 * Scales a vec2 by a scalar number
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the vector to scale
 * @param {Number} b amount to scale the vector by
 * @returns {vec2} out
 */
vec2_p vec2_scale(vec2_p out, vec2_cp a, scalar_t b) {
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	return out;
};

/**
 * Adds two vec2's after scaling the second operand by a scalar value
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @param {Number} scale the amount to scale b by before adding
 * @returns {vec2} out
 */
vec2_p vec2_madd(vec2_p out, vec2_cp a, vec2_cp b, scalar_t scale) {
	out[0] = a[0] + (b[0] * scale);
	out[1] = a[1] + (b[1] * scale);
	return out;
};

/**
 * Calculates the euclidian distance between two vec2's
 *
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {Number} distance between a and b
 */
scalar_t vec2_dist(vec2_cp a, vec2_cp b) {
	scalar_t x = b[0] - a[0], y = b[1] - a[1];
	return sqrt(x*x + y*y);
};

/**
 * Calculates the squared euclidian distance between two vec2's
 *
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {Number} squared distance between a and b
 */
scalar_t vec2_dist2(vec2_cp a, vec2_cp b) {
	scalar_t x = b[0] - a[0], y = b[1] - a[1];
	return x*x + y*y;
};

/**
 * Calculates the length of a vec2
 *
 * @param {vec2} a vector to calculate length of
 * @returns {Number} length of a
 */
scalar_t vec2_len(vec2_cp a) {
	scalar_t x = a[0], y = a[1];
	return sqrt(x*x + y*y);
};

/**
 * Calculates the squared length of a vec2
 *
 * @param {vec2} a vector to calculate squared length of
 * @returns {Number} squared length of a
 */
scalar_t vec2_len2(vec2_cp a) {
	scalar_t x = a[0], y = a[1];
	return x*x + y*y;
};

/**
 * Negates the components of a vec2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a vector to negate
 * @returns {vec2} out
 */
vec2_p vec2_negate(vec2_p out, vec2_cp a) {
	out[0] = -a[0];
	out[1] = -a[1];
	return out;
};

/**
 * Returns the inverse of the components of a vec2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a vector to invert
 * @returns {vec2} out
 */
vec2_p vec2_inverse(vec2_p out, vec2_cp a) {
	out[0] = 1.0 / a[0];
	out[1] = 1.0 / a[1];
	return out;
};

/**
 * Normalize a vec2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a vector to normalize
 * @returns {vec2} out
 */
vec2_p normalize(vec2_p out, vec2_cp a) {
	scalar_t x = a[0], y = a[1];
	scalar_t len = x*x + y*y;
	if (len > 0) {
		//TODO: evaluate use of glm_invsqrt here?
		len = 1 / sqrt(len);
		out[0] = a[0] * len;
		out[1] = a[1] * len;
	}
	return out;
};

/**
 * Calculates the dot product of two vec2's
 *
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {Number} dot product of a and b
 */
scalar_t vec2_dot(vec2_cp a, vec2_cp b) {
	return a[0] * b[0] + a[1] * b[1];
};

/**
 * Computes the cross product of two vec2's
 * Note that the cross product must by definition produce a 3D vector
 *
 * @param {vec3} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @returns {vec3} out
 */
vec3_p vec2_cross(vec3_p out, vec2_cp a, vec2_cp b) {
	scalar_t z = a[0] * b[1] - a[1] * b[0];
	out[0] = out[1] = 0;
	out[2] = z;
	return out;
};

/**
 * Performs a linear interpolation between two vec2's
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the first operand
 * @param {vec2} b the second operand
 * @param {Number} t interpolation amount between the two inputs
 * @returns {vec2} out
 */
vec2_p vec2_lerp(vec2_p out, vec2_cp a, vec2_cp b, scalar_t t) {
	scalar_t ax = a[0], ay = a[1];
	out[0] = ax + t * (b[0] - ax);
	out[1] = ay + t * (b[1] - ay);
	return out;
};

/**
 * Generates a random vector with the given scale
 *
 * @param {vec2} out the receiving vector
 * @param {Number} [scale] Length of the resulting vector. If ommitted, a unit vector will be returned
 * @returns {vec2} out
 */
/*
export function random(out, scale) {
	scale = scale || 1.0;
	var r = glMatrix.RANDOM() * 2.0 * Math.PI;
	out[0] = Math.cos(r) * scale;
	out[1] = Math.sin(r) * scale;
	return out;
};
*/

/**
 * Transforms the vec2 with a mat2
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the vector to transform
 * @param {mat2} m matrix to transform with
 * @returns {vec2} out
 */
vec2_p vec2_xform_mat2(vec2_p out, vec2_cp a, mat2_cp m) {
	scalar_t x = a[0], y = a[1];
	out[0] = m[0] * x + m[2] * y;
	out[1] = m[1] * x + m[3] * y;
	return out;
};

/**
 * Transforms the vec2 with a mat2d
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the vector to transform
 * @param {mat2d} m matrix to transform with
 * @returns {vec2} out
 */
vec2_p vec2_xform_mat2d(vec2_p out, vec2_cp a, mat2d_cp m) {
	scalar_t x = a[0], y = a[1];
	out[0] = m[0] * x + m[2] * y + m[4];
	out[1] = m[1] * x + m[3] * y + m[5];
	return out;
};

/**
 * Transforms the vec2 with a mat3
 * 3rd vector component is implicitly '1'
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the vector to transform
 * @param {mat3} m matrix to transform with
 * @returns {vec2} out
 */
vec2_p vec2_xform_mat3(vec2_p out, vec2_cp a, mat3_cp m) {
	scalar_t x = a[0], y = a[1];
	out[0] = m[0] * x + m[3] * y + m[6];
	out[1] = m[1] * x + m[4] * y + m[7];
	return out;
};

/**
 * Transforms the vec2 with a mat4
 * 3rd vector component is implicitly '0'
 * 4th vector component is implicitly '1'
 *
 * @param {vec2} out the receiving vector
 * @param {vec2} a the vector to transform
 * @param {mat4} m matrix to transform with
 * @returns {vec2} out
 */
vec2_p vec2_xform_mat4(vec2_p out, vec2_cp a, mat4_cp m) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	out[0] = m[0] * x + m[4] * y + m[12];
	out[1] = m[1] * x + m[5] * y + m[13];
	return out;
}

/**
 * Returns whether or not the vectors exactly have the same elements in the same position (when compared with ===)
 *
 * @param {vec2} a The first vector.
 * @param {vec2} b The second vector.
 * @returns {Boolean} True if the vectors are equal, false otherwise.
 */
int vec2_exact_equals(vec2_cp a, vec2_cp b) {
	return a[0] == b[0] && a[1] == b[1];
}

/**
 * Returns whether or not the vectors have approximately the same elements in the same position.
 *
 * @param {vec2} a The first vector.
 * @param {vec2} b The second vector.
 * @returns {Boolean} True if the vectors are equal, false otherwise.
 */
int vec2_equals(vec2_cp a, vec2_cp b) {
	scalar_t a0 = a[0], a1 = a[1];
	scalar_t b0 = b[0], b1 = b[1];
	return (fabs(a0 - b0) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a0), fabs(b0))) &&
		fabs(a1 - b1) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a1), fabs(b1))));
}
