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

/* C port by Traian Stanev */

#include "defs.h"

/**
 * 4 Dimensional Vector
 * @module vec4
 */

/**
 * Copy the values from one vec4 to another
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the source vector
 * @returns {vec4} out
 */
vec4_p vec4_copy(vec4_p out, vec4_cp a) {
	out[0] = a[0];
	out[1] = a[1];
	out[2] = a[2];
	out[3] = a[3];
	return out;
}

/**
 * Set the components of a vec4 to the given values
 *
 * @param {vec4} out the receiving vector
 * @param {Number} x X component
 * @param {Number} y Y component
 * @param {Number} z Z component
 * @param {Number} w W component
 * @returns {vec4} out
 */
vec4_p vec4_set(vec4_p out, scalar_t x, scalar_t y, scalar_t z, scalar_t w) {
	out[0] = x;
	out[1] = y;
	out[2] = z;
	out[3] = w;
	return out;
}

/**
 * Adds two vec4's
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {vec4} out
 */
vec4_p vec4_add(vec4_p out, vec4_cp a, vec4_cp b) {
	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
	out[2] = a[2] + b[2];
	out[3] = a[3] + b[3];
	return out;
}

/**
 * Subtracts vector b from vector a
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {vec4} out
 */
vec4_p vec4_sub(vec4_p out, vec4_cp a, vec4_cp b) {
	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
	out[2] = a[2] - b[2];
	out[3] = a[3] - b[3];
	return out;
}

/**
 * Multiplies two vec4's
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {vec4} out
 */
vec4_p vec4_mul(vec4_p out, vec4_cp a, vec4_cp b) {
	out[0] = a[0] * b[0];
	out[1] = a[1] * b[1];
	out[2] = a[2] * b[2];
	out[3] = a[3] * b[3];
	return out;
}

/**
 * Divides two vec4's
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {vec4} out
 */
vec4_p vec4_div(vec4_p out, vec4_cp a, vec4_cp b) {
	out[0] = a[0] / b[0];
	out[1] = a[1] / b[1];
	out[2] = a[2] / b[2];
	out[3] = a[3] / b[3];
	return out;
}

/**
 * Math.ceil the components of a vec4
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a vector to ceil
 * @returns {vec4} out
 */
vec4_p vec4_ceil(vec4_p out, vec4_cp a) {
	out[0] = ceil(a[0]);
	out[1] = ceil(a[1]);
	out[2] = ceil(a[2]);
	out[3] = ceil(a[3]);
	return out;
}

/**
 * Math.floor the components of a vec4
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a vector to floor
 * @returns {vec4} out
 */
vec4_p vec4_floor(vec4_p out, vec4_cp a) {
	out[0] = floor(a[0]);
	out[1] = floor(a[1]);
	out[2] = floor(a[2]);
	out[3] = floor(a[3]);
	return out;
}

/**
 * Returns the minimum of two vec4's
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {vec4} out
 */
vec4_p vec4_min(vec4_p out, vec4_cp a, vec4_cp b) {
	out[0] = fmin(a[0], b[0]);
	out[1] = fmin(a[1], b[1]);
	out[2] = fmin(a[2], b[2]);
	out[3] = fmin(a[3], b[3]);
	return out;
}

/**
 * Returns the maximum of two vec4's
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {vec4} out
 */
vec4_p vec4_max(vec4_p out, vec4_cp a, vec4_cp b) {
	out[0] = fmax(a[0], b[0]);
	out[1] = fmax(a[1], b[1]);
	out[2] = fmax(a[2], b[2]);
	out[3] = fmax(a[3], b[3]);
	return out;
}

/**
 * Math.round the components of a vec4
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a vector to round
 * @returns {vec4} out
 */
vec4_p vec4_round(vec4_p out, vec4_cp a) {
	out[0] = round(a[0]);
	out[1] = round(a[1]);
	out[2] = round(a[2]);
	out[3] = round(a[3]);
	return out;
}

/**
 * Scales a vec4 by a scalar number
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the vector to scale
 * @param {Number} b amount to scale the vector by
 * @returns {vec4} out
 */
vec4_p vec4_scale(vec4_p out, vec4_cp a, scalar_t b) {
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;
	out[3] = a[3] * b;
	return out;
}

/**
 * Adds two vec4's after scaling the second operand by a scalar value
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @param {Number} scale the amount to scale b by before adding
 * @returns {vec4} out
 */
vec4_p vec4_madd(vec4_p out, vec4_cp a, vec4_cp b, scalar_t scale) {
	out[0] = a[0] + (b[0] * scale);
	out[1] = a[1] + (b[1] * scale);
	out[2] = a[2] + (b[2] * scale);
	out[3] = a[3] + (b[3] * scale);
	return out;
}

/**
 * Calculates the euclidian distance between two vec4's
 *
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {Number} distance between a and b
 */
scalar_t vec4_dist(vec4_cp a, vec4_cp b) {
	scalar_t x = b[0] - a[0];
	scalar_t y = b[1] - a[1];
	scalar_t z = b[2] - a[2];
	scalar_t w = b[3] - a[3];
	return sqrt(x*x + y*y + z*z + w*w);
}

/**
 * Calculates the squared euclidian distance between two vec4's
 *
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {Number} squared distance between a and b
 */
scalar_t vec4_dist2(vec4_cp a, vec4_cp b) {
	scalar_t x = b[0] - a[0];
	scalar_t y = b[1] - a[1];
	scalar_t z = b[2] - a[2];
	scalar_t w = b[3] - a[3];
	return x*x + y*y + z*z + w*w;
}

/**
 * Calculates the length of a vec4
 *
 * @param {vec4} a vector to calculate length of
 * @returns {Number} length of a
 */
scalar_t vec4_len(vec4_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	scalar_t w = a[3];
	return sqrt(x*x + y*y + z*z + w*w);
}

/**
 * Calculates the squared length of a vec4
 *
 * @param {vec4} a vector to calculate squared length of
 * @returns {Number} squared length of a
 */
scalar_t vec4_len2(vec4_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	scalar_t w = a[3];
	return x*x + y*y + z*z + w*w;
}

/**
 * Negates the components of a vec4
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a vector to negate
 * @returns {vec4} out
 */
vec4_p vec4_negate(vec4_p out, vec4_cp a) {
	out[0] = -a[0];
	out[1] = -a[1];
	out[2] = -a[2];
	out[3] = -a[3];
	return out;
}

/**
 * Returns the inverse of the components of a vec4
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a vector to invert
 * @returns {vec4} out
 */
vec4_p vec4_inverse(vec4_p out, vec4_cp a) {
	out[0] = 1.0 / a[0];
	out[1] = 1.0 / a[1];
	out[2] = 1.0 / a[2];
	out[3] = 1.0 / a[3];
	return out;
}

/**
 * Normalize a vec4
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a vector to normalize
 * @returns {vec4} out
 */
vec4_p vec4_normalize(vec4_p out, vec4_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	scalar_t w = a[3];
	scalar_t len = x*x + y*y + z*z + w*w;
	if (len > 0) {
		len = 1 / sqrt(len);
		out[0] = x * len;
		out[1] = y * len;
		out[2] = z * len;
		out[3] = w * len;
	}
	return out;
}

/**
 * Calculates the dot product of two vec4's
 *
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @returns {Number} dot product of a and b
 */
scalar_t vec4_dot(vec4_cp a, vec4_cp b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

/**
 * Performs a linear interpolation between two vec4's
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the first operand
 * @param {vec4} b the second operand
 * @param {Number} t interpolation amount between the two inputs
 * @returns {vec4} out
 */
vec4_p vec4_lerp(vec4_p out, vec4_cp a, vec4_cp b, scalar_t t) {
	scalar_t ax = a[0];
	scalar_t ay = a[1];
	scalar_t az = a[2];
	scalar_t aw = a[3];
	out[0] = ax + t * (b[0] - ax);
	out[1] = ay + t * (b[1] - ay);
	out[2] = az + t * (b[2] - az);
	out[3] = aw + t * (b[3] - aw);
	return out;
}

/**
 * Generates a random vector with the given scale
 *
 * @param {vec4} out the receiving vector
 * @param {Number} [scale] Length of the resulting vector. If ommitted, a unit vector will be returned
 * @returns {vec4} out
 */
/*
export function random(out, vectorScale) {
	vectorScale = vectorScale || 1.0;

	//TODO: This is a pretty awful way of doing this. Find something better.
	out[0] = glMatrix.RANDOM();
	out[1] = glMatrix.RANDOM();
	out[2] = glMatrix.RANDOM();
	out[3] = glMatrix.RANDOM();
	normalize(out, out);
	scale(out, out, vectorScale);
	return out;
}
*/

/**
 * Transforms the vec4 with a mat4.
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the vector to transform
 * @param {mat4} m matrix to transform with
 * @returns {vec4} out
 */
vec4_p vec4_xform_mat4(vec4_p out, vec4_cp a, mat4_cp m) {
	scalar_t x = a[0], y = a[1], z = a[2], w = a[3];
	out[0] = m[0] * x + m[4] * y + m[8] * z + m[12] * w;
	out[1] = m[1] * x + m[5] * y + m[9] * z + m[13] * w;
	out[2] = m[2] * x + m[6] * y + m[10] * z + m[14] * w;
	out[3] = m[3] * x + m[7] * y + m[11] * z + m[15] * w;
	return out;
}

/**
 * Transforms the vec4 with a quat
 *
 * @param {vec4} out the receiving vector
 * @param {vec4} a the vector to transform
 * @param {quat} q quaternion to transform with
 * @returns {vec4} out
 */
vec4_p vec4_xform_quat(vec4_p out, vec4_cp a, quat_cp q) {
	scalar_t x = a[0], y = a[1], z = a[2];
	scalar_t qx = q[0], qy = q[1], qz = q[2], qw = q[3];

	// calculate quat * vec
	scalar_t ix = qw * x + qy * z - qz * y;
	scalar_t iy = qw * y + qz * x - qx * z;
	scalar_t iz = qw * z + qx * y - qy * x;
	scalar_t iw = -qx * x - qy * y - qz * z;

	// calculate result * inverse quat
	out[0] = ix * qw + iw * -qx + iy * -qz - iz * -qy;
	out[1] = iy * qw + iw * -qy + iz * -qx - ix * -qz;
	out[2] = iz * qw + iw * -qz + ix * -qy - iy * -qx;
	out[3] = a[3];
	return out;
}

/**
 * Returns whether or not the vectors have exactly the same elements in the same position (when compared with ===)
 *
 * @param {vec4} a The first vector.
 * @param {vec4} b The second vector.
 * @returns {Boolean} True if the vectors are equal, false otherwise.
 */
int vec4_exact_equals(vec4_cp a, vec4_cp b) {
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

/**
 * Returns whether or not the vectors have approximately the same elements in the same position.
 *
 * @param {vec4} a The first vector.
 * @param {vec4} b The second vector.
 * @returns {Boolean} True if the vectors are equal, false otherwise.
 */
int vec4_equals(vec4_cp a, vec4_cp b) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
	scalar_t b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
	return (fabs(a0 - b0) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a0), fabs(b0))) &&
		fabs(a1 - b1) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a1), fabs(b1))) &&
		fabs(a2 - b2) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a2), fabs(b2))) &&
		fabs(a3 - b3) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a3), fabs(b3))));
}
