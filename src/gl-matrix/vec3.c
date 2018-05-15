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

#pragma once 

#include "defs.h"

/**
 * 3 Dimensional Vector
 * @module vec3
 */

/**
 * Calculates the length of a vec3
 *
 * @param {vec3} a vector to calculate length of
 * @returns {Number} length of a
 */
scalar_t vec3_length(vec3_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	return sqrt(x*x + y*y + z*z);
}

/**
 * Copy the values from one vec3 to another
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the source vector
 * @returns {vec3} out
 */
vec3_p vec3_copy(vec3_p out, vec3_cp a) {
	out[0] = a[0];
	out[1] = a[1];
	out[2] = a[2];
	return out;
}

/**
 * Set the components of a vec3 to the given values
 *
 * @param {vec3} out the receiving vector
 * @param {Number} x X component
 * @param {Number} y Y component
 * @param {Number} z Z component
 * @returns {vec3} out
 */
vec3_p vec3_set(vec3_p out, scalar_t x, scalar_t y, scalar_t z) {
	out[0] = x;
	out[1] = y;
	out[2] = z;
	return out;
}

/**
 * Adds two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_add(vec3_p out, vec3_cp a, vec3_cp b) {
	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
	out[2] = a[2] + b[2];
	return out;
}

/**
 * Subtracts vector b from vector a
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_sub(vec3_p out, vec3_cp a, vec3_cp b) {
	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
	out[2] = a[2] - b[2];
	return out;
}

/**
 * Multiplies two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_mul(vec3_p out, vec3_cp a, vec3_cp b) {
	out[0] = a[0] * b[0];
	out[1] = a[1] * b[1];
	out[2] = a[2] * b[2];
	return out;
}

/**
 * Divides two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_div(vec3_p out, vec3_cp a, vec3_cp b) {
	out[0] = a[0] / b[0];
	out[1] = a[1] / b[1];
	out[2] = a[2] / b[2];
	return out;
}

/**
 * Math.ceil the components of a vec3
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a vector to ceil
 * @returns {vec3} out
 */
vec3_p vec3_ceil(vec3_p out, vec3_cp a) {
	out[0] = ceil(a[0]);
	out[1] = ceil(a[1]);
	out[2] = ceil(a[2]);
	return out;
}

/**
 * Math.floor the components of a vec3
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a vector to floor
 * @returns {vec3} out
 */
vec3_p vec3_floor(vec3_p out, vec3_cp a) {
	out[0] = floor(a[0]);
	out[1] = floor(a[1]);
	out[2] = floor(a[2]);
	return out;
}

/**
 * Returns the minimum of two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_min(vec3_p out, vec3_cp a, vec3_cp b) {
	out[0] = fmin(a[0], b[0]);
	out[1] = fmin(a[1], b[1]);
	out[2] = fmin(a[2], b[2]);
	return out;
}

/**
 * Returns the maximum of two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_max(vec3_p out, vec3_cp a, vec3_cp b) {
	out[0] = fmax(a[0], b[0]);
	out[1] = fmax(a[1], b[1]);
	out[2] = fmax(a[2], b[2]);
	return out;
}

/**
 * Math.round the components of a vec3
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a vector to round
 * @returns {vec3} out
 */
vec3_p vec3_round(vec3_p out, vec3_cp a) {
	out[0] = round(a[0]);
	out[1] = round(a[1]);
	out[2] = round(a[2]);
	return out;
}

/**
 * Scales a vec3 by a scalar number
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the vector to scale
 * @param {Number} b amount to scale the vector by
 * @returns {vec3} out
 */
vec3_p vec3_scale(vec3_p out, vec3_cp a, scalar_t b) {
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;
	return out;
}

/**
 * Adds two vec3's after scaling the second operand by a scalar value
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @param {Number} scale the amount to scale b by before adding
 * @returns {vec3} out
 */
vec3_p vec3_madd(vec3_p out, vec3_cp a, vec3_cp b, scalar_t scale) {
	out[0] = a[0] + (b[0] * scale);
	out[1] = a[1] + (b[1] * scale);
	out[2] = a[2] + (b[2] * scale);
	return out;
}

/**
 * Calculates the euclidian distance between two vec3's
 *
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {Number} distance between a and b
 */
scalar_t vec3_dist(vec3_cp a, vec3_cp b) {
	scalar_t x = b[0] - a[0];
	scalar_t y = b[1] - a[1];
	scalar_t z = b[2] - a[2];
	return sqrt(x*x + y*y + z*z);
}

/**
 * Calculates the squared euclidian distance between two vec3's
 *
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {Number} squared distance between a and b
 */
scalar_t vec3_dist2(vec3_cp a, vec3_cp b) {
	scalar_t x = b[0] - a[0];
	scalar_t y = b[1] - a[1];
	scalar_t z = b[2] - a[2];
	return x*x + y*y + z*z;
}

/**
 * Calculates the length of a vec3
 *
 * @param {vec3} a vector to calculate length of
 * @returns {Number} squared length of a
 */
scalar_t vec3_len(vec3_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	return sqrt(x*x + y*y + z*z);
}

/**
 * Calculates the squared length of a vec3
 *
 * @param {vec3} a vector to calculate squared length of
 * @returns {Number} squared length of a
 */
scalar_t vec3_len2(vec3_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	return x*x + y*y + z*z;
}

/**
 * Negates the components of a vec3
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a vector to negate
 * @returns {vec3} out
 */
vec3_p vec3_negate(vec3_p out, vec3_cp a) {
	out[0] = -a[0];
	out[1] = -a[1];
	out[2] = -a[2];
	return out;
}

/**
 * Returns the inverse of the components of a vec3
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a vector to invert
 * @returns {vec3} out
 */
vec3_p vec3_inverse(vec3_p out, vec3_cp a) {
	out[0] = 1.0 / a[0];
	out[1] = 1.0 / a[1];
	out[2] = 1.0 / a[2];
	return out;
}

/**
 * Normalize a vec3
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a vector to normalize
 * @returns {vec3} out
 */
vec3_p vec3_normalize(vec3_p out, vec3_cp a) {
	scalar_t x = a[0];
	scalar_t y = a[1];
	scalar_t z = a[2];
	scalar_t len = x*x + y*y + z*z;
	if (len > 0) {
		//TODO: evaluate use of glm_invsqrt here?
		len = 1 / sqrt(len);
		out[0] = a[0] * len;
		out[1] = a[1] * len;
		out[2] = a[2] * len;
	}
	return out;
}

/**
 * Calculates the dot product of two vec3's
 *
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {Number} dot product of a and b
 */
scalar_t vec3_dot(vec3_cp a, vec3_cp b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Computes the cross product of two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @returns {vec3} out
 */
vec3_p vec3_cross(vec3_p out, vec3_cp a, vec3_cp b) {
	scalar_t ax = a[0], ay = a[1], az = a[2];
	scalar_t bx = b[0], by = b[1], bz = b[2];

	out[0] = ay * bz - az * by;
	out[1] = az * bx - ax * bz;
	out[2] = ax * by - ay * bx;
	return out;
}

/**
 * Performs a linear interpolation between two vec3's
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @param {Number} t interpolation amount between the two inputs
 * @returns {vec3} out
 */
vec3_p vec3_lerp(vec3_p out, vec3_cp a, vec3_cp b, scalar_t t) {
	scalar_t ax = a[0];
	scalar_t ay = a[1];
	scalar_t az = a[2];
	out[0] = ax + t * (b[0] - ax);
	out[1] = ay + t * (b[1] - ay);
	out[2] = az + t * (b[2] - az);
	return out;
}

/**
 * Performs a hermite interpolation with two control points
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @param {vec3} c the third operand
 * @param {vec3} d the fourth operand
 * @param {Number} t interpolation amount between the two inputs
 * @returns {vec3} out
 */
vec3_p vec3_hermite(vec3_p out, vec3_cp a, vec3_cp b, vec3_cp c, vec3_cp d, scalar_t t) {
	scalar_t factorTimes2 = t * t;
	scalar_t factor1 = factorTimes2 * (2 * t - 3) + 1;
	scalar_t factor2 = factorTimes2 * (t - 2) + t;
	scalar_t factor3 = factorTimes2 * (t - 1);
	scalar_t factor4 = factorTimes2 * (3 - 2 * t);

	out[0] = a[0] * factor1 + b[0] * factor2 + c[0] * factor3 + d[0] * factor4;
	out[1] = a[1] * factor1 + b[1] * factor2 + c[1] * factor3 + d[1] * factor4;
	out[2] = a[2] * factor1 + b[2] * factor2 + c[2] * factor3 + d[2] * factor4;

	return out;
}

/**
 * Performs a bezier interpolation with two control points
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the first operand
 * @param {vec3} b the second operand
 * @param {vec3} c the third operand
 * @param {vec3} d the fourth operand
 * @param {Number} t interpolation amount between the two inputs
 * @returns {vec3} out
 */
vec3_p vec3_bezier(vec3_p out, vec3_cp a, vec3_cp b, vec3_cp c, vec3_cp d, scalar_t t) {
	scalar_t inverseFactor = 1 - t;
	scalar_t inverseFactorTimesTwo = inverseFactor * inverseFactor;
	scalar_t factorTimes2 = t * t;
	scalar_t factor1 = inverseFactorTimesTwo * inverseFactor;
	scalar_t factor2 = 3 * t * inverseFactorTimesTwo;
	scalar_t factor3 = 3 * factorTimes2 * inverseFactor;
	scalar_t factor4 = factorTimes2 * t;

	out[0] = a[0] * factor1 + b[0] * factor2 + c[0] * factor3 + d[0] * factor4;
	out[1] = a[1] * factor1 + b[1] * factor2 + c[1] * factor3 + d[1] * factor4;
	out[2] = a[2] * factor1 + b[2] * factor2 + c[2] * factor3 + d[2] * factor4;

	return out;
}

/**
 * Generates a random vector with the given scale
 *
 * @param {vec3} out the receiving vector
 * @param {Number} [scale] Length of the resulting vector. If ommitted, a unit vector will be returned
 * @returns {vec3} out
 */
/*
export function random(out, scale) {
	scale = scale || 1.0;

	let r = glMatrix.RANDOM() * 2.0 * Math.PI;
	let z = (glMatrix.RANDOM() * 2.0) - 1.0;
	let zScale = Math.sqrt(1.0-z*z) * scale;

	out[0] = Math.cos(r) * zScale;
	out[1] = Math.sin(r) * zScale;
	out[2] = z * scale;
	return out;
}
*/

/**
 * Transforms the vec3 with a mat4.
 * 4th vector component is implicitly '1'
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the vector to transform
 * @param {mat4} m matrix to transform with
 * @returns {vec3} out
 */
vec3_p vec3_xform_mat4(vec3_p out, vec3_cp a, mat4_cp m) {
	scalar_t x = a[0], y = a[1], z = a[2];
	scalar_t w = m[3] * x + m[7] * y + m[11] * z + m[15];
	if (w == 0) {
		w = 1;
	}
	out[0] = (m[0] * x + m[4] * y + m[8] * z + m[12]) / w;
	out[1] = (m[1] * x + m[5] * y + m[9] * z + m[13]) / w;
	out[2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) / w;
	return out;
}

/**
 * Transforms the vec3 with a mat3.
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the vector to transform
 * @param {mat3} m the 3x3 matrix to transform with
 * @returns {vec3} out
 */
vec3_p vec3_xform_mat3(vec3_p out, vec3_cp a, mat3_cp m) {
	scalar_t x = a[0], y = a[1], z = a[2];
	out[0] = x * m[0] + y * m[3] + z * m[6];
	out[1] = x * m[1] + y * m[4] + z * m[7];
	out[2] = x * m[2] + y * m[5] + z * m[8];
	return out;
}

/**
 * Transforms the vec3 with a quat
 * Can also be used for dual quaternions. (Multiply it with the real part)
 *
 * @param {vec3} out the receiving vector
 * @param {vec3} a the vector to transform
 * @param {quat} q quaternion to transform with
 * @returns {vec3} out
 */
vec3_p vec3_xform_quat(vec3_p out, vec3_cp a, quat_cp q) {
		// benchmarks: https://jsperf.com/quaternion-transform-vec3-implementations-fixed
		scalar_t qx = q[0], qy = q[1], qz = q[2], qw = q[3];
		scalar_t x = a[0], y = a[1], z = a[2];
		// var qvec = [qx, qy, qz];
		// var uv = vec3.cross([], qvec, a);
		scalar_t uvx = qy * z - qz * y,
			uvy = qz * x - qx * z,
			uvz = qx * y - qy * x;
		// var uuv = vec3.cross([], qvec, uv);
		scalar_t uuvx = qy * uvz - qz * uvy,
			uuvy = qz * uvx - qx * uvz,
			uuvz = qx * uvy - qy * uvx;
		// vec3.scale(uv, uv, 2 * w);
		scalar_t w2 = qw * 2;
		uvx *= w2;
		uvy *= w2;
		uvz *= w2;
		// vec3.scale(uuv, uuv, 2);
		uuvx *= 2;
		uuvy *= 2;
		uuvz *= 2;
		// return vec3.add(out, a, vec3.add(out, uv, uuv));
		out[0] = x + uvx + uuvx;
		out[1] = y + uvy + uuvy;
		out[2] = z + uvz + uuvz;
		return out;
}

/**
 * Rotate a 3D vector around the x-axis
 * @param {vec3} out The receiving vec3
 * @param {vec3} a The vec3 point to rotate
 * @param {vec3} b The origin of the rotation
 * @param {Number} c The angle of rotation
 * @returns {vec3} out
 */
vec3_p vec3_rotate_x(vec3_p out, vec3_cp a, vec3_cp b, scalar_t c){
	scalar_t p[3], r[3];
	//Translate point to the origin
	p[0] = a[0] - b[0];
	p[1] = a[1] - b[1];
	p[2] = a[2] - b[2];

	//perform rotation
	r[0] = p[0];
	r[1] = p[1]*cos(c) - p[2]*sin(c);
	r[2] = p[1]*sin(c) + p[2]*cos(c);

	//translate to correct position
	out[0] = r[0] + b[0];
	out[1] = r[1] + b[1];
	out[2] = r[2] + b[2];

	return out;
}

/**
 * Rotate a 3D vector around the y-axis
 * @param {vec3} out The receiving vec3
 * @param {vec3} a The vec3 point to rotate
 * @param {vec3} b The origin of the rotation
 * @param {Number} c The angle of rotation
 * @returns {vec3} out
 */
vec3_p vec3_rotate_y(vec3_p out, vec3_cp a, vec3_cp b, scalar_t c){
	scalar_t p[3], r[3];
	//Translate point to the origin
	p[0] = a[0] - b[0];
	p[1] = a[1] - b[1];
	p[2] = a[2] - b[2];

	//perform rotation
	r[0] = p[2]*sin(c) + p[0]*cos(c);
	r[1] = p[1];
	r[2] = p[2]*cos(c) - p[0]*sin(c);

	//translate to correct position
	out[0] = r[0] + b[0];
	out[1] = r[1] + b[1];
	out[2] = r[2] + b[2];

	return out;
}

/**
 * Rotate a 3D vector around the z-axis
 * @param {vec3} out The receiving vec3
 * @param {vec3} a The vec3 point to rotate
 * @param {vec3} b The origin of the rotation
 * @param {Number} c The angle of rotation
 * @returns {vec3} out
 */
vec3_p vec3_rotate_z(vec3_p out, vec3_cp a, vec3_cp b, scalar_t c){
	scalar_t p[3], r[3];
	//Translate point to the origin
	p[0] = a[0] - b[0];
	p[1] = a[1] - b[1];
	p[2] = a[2] - b[2];

	//perform rotation
	r[0] = p[0]*cos(c) - p[1]*sin(c);
	r[1] = p[0]*sin(c) + p[1]*cos(c);
	r[2] = p[2];

	//translate to correct position
	out[0] = r[0] + b[0];
	out[1] = r[1] + b[1];
	out[2] = r[2] + b[2];

	return out;
}

/**
 * Get the angle between two 3D vectors
 * @param {vec3} a The first operand
 * @param {vec3} b The second operand
 * @returns {Number} The angle in radians
 */
scalar_t vec3_angle(vec3_cp a, vec3_cp b) {
	scalar_t tempA[] = {a[0], a[1], a[2]};
	scalar_t tempB[] = {b[0], b[1], b[2]};

	vec3_normalize(tempA, tempA);
	vec3_normalize(tempB, tempB);

	scalar_t cosine = vec3_dot(tempA, tempB);

	if(cosine > 1.0) {
		return 0;
	}
	else if(cosine < -1.0) {
		return GLM_PI;
	} else {
		return acos(cosine);
	}
}

/**
 * Returns whether or not the vectors have exactly the same elements in the same position (when compared with ===)
 *
 * @param {vec3} a The first vector.
 * @param {vec3} b The second vector.
 * @returns {Boolean} True if the vectors are equal, false otherwise.
 */
int vec3_exact_equals(vec3_cp a, vec3_cp b) {
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

/**
 * Returns whether or not the vectors have approximately the same elements in the same position.
 *
 * @param {vec3} a The first vector.
 * @param {vec3} b The second vector.
 * @returns {Boolean} True if the vectors are equal, false otherwise.
 */
int vec3_equals(vec3_cp a, vec3_cp b) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2];
	scalar_t b0 = b[0], b1 = b[1], b2 = b[2];
	return (fabs(a0 - b0) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a0), fabs(b0))) &&
		fabs(a1 - b1) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a1), fabs(b1))) &&
		fabs(a2 - b2) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a2), fabs(b2))));
}
