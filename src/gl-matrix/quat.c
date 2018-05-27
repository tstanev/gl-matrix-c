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
#include "vec3.c"
#include "vec4.c"

/**
 * Quaternion
 * @module quat
 */

/**
 * Set a quat to the identity quaternion
 *
 * @param {quat} out the receiving quaternion
 * @returns {quat} out
 */
quat_p quat_identity(quat_p out) {
	out[0] = 0;
	out[1] = 0;
	out[2] = 0;
	out[3] = 1;
	return out;
}

/**
 * Sets a quat from the given angle and rotation axis,
 * then returns it.
 *
 * @param {quat} out the receiving quaternion
 * @param {vec3} axis the axis around which to rotate
 * @param {Number} rad the angle in radians
 * @returns {quat} out
 **/
quat_p quat_set_axis_angle(quat_p out, vec3_cp axis, scalar_t rad) {
	rad = rad * 0.5;
	scalar_t s = sin(rad);
	out[0] = s * axis[0];
	out[1] = s * axis[1];
	out[2] = s * axis[2];
	out[3] = cos(rad);
	return out;
}

/**
 * Gets the rotation axis and angle for a given
 *  quaternion. If a quaternion is created with
 *  setAxisAngle, this method will return the same
 *  values as providied in the original parameter list
 *  OR functionally equivalent values.
 * Example: The quaternion formed by axis [0, 0, 1] and
 *  angle -90 is the same as the quaternion formed by
 *  [0, 0, 1] and 270. This method favors the latter.
 * @param  {vec3} out_axis  Vector receiving the axis of rotation
 * @param  {quat} q     Quaternion to be decomposed
 * @return {Number}     Angle, in radians, of the rotation
 */
scalar_t quat_get_axis_angle(vec3_p out_axis, quat_cp q) {
	scalar_t rad = acos(q[3]) * 2.0;
	scalar_t s = sin(rad / 2.0);
	if (s != 0.0) {
		out_axis[0] = q[0] / s;
		out_axis[1] = q[1] / s;
		out_axis[2] = q[2] / s;
	} else {
		// If s is zero, return any axis (no rotation - axis does not matter)
		out_axis[0] = 1;
		out_axis[1] = 0;
		out_axis[2] = 0;
	}
	return rad;
}

/**
 * Multiplies two quat's
 *
 * @param {quat} out the receiving quaternion
 * @param {quat} a the first operand
 * @param {quat} b the second operand
 * @returns {quat} out
 */
quat_p quat_multiply(quat_p out, quat_cp a, quat_cp b) {
	scalar_t ax = a[0], ay = a[1], az = a[2], aw = a[3];
	scalar_t bx = b[0], by = b[1], bz = b[2], bw = b[3];

	out[0] = ax * bw + aw * bx + ay * bz - az * by;
	out[1] = ay * bw + aw * by + az * bx - ax * bz;
	out[2] = az * bw + aw * bz + ax * by - ay * bx;
	out[3] = aw * bw - ax * bx - ay * by - az * bz;
	return out;
}

/**
 * Rotates a quaternion by the given angle about the X axis
 *
 * @param {quat} out quat receiving operation result
 * @param {quat} a quat to rotate
 * @param {number} rad angle (in radians) to rotate
 * @returns {quat} out
 */
quat_p quat_rotate_x(quat_p out, quat_cp a, scalar_t rad) {
	rad *= 0.5;

	scalar_t ax = a[0], ay = a[1], az = a[2], aw = a[3];
	scalar_t bx = sin(rad), bw = cos(rad);

	out[0] = ax * bw + aw * bx;
	out[1] = ay * bw + az * bx;
	out[2] = az * bw - ay * bx;
	out[3] = aw * bw - ax * bx;
	return out;
}

/**
 * Rotates a quaternion by the given angle about the Y axis
 *
 * @param {quat} out quat receiving operation result
 * @param {quat} a quat to rotate
 * @param {number} rad angle (in radians) to rotate
 * @returns {quat} out
 */
quat_p quat_rotate_y(quat_p out, quat_cp a, scalar_t rad) {
	rad *= 0.5;

	scalar_t ax = a[0], ay = a[1], az = a[2], aw = a[3];
	scalar_t by = sin(rad), bw = cos(rad);

	out[0] = ax * bw - az * by;
	out[1] = ay * bw + aw * by;
	out[2] = az * bw + ax * by;
	out[3] = aw * bw - ay * by;
	return out;
}

/**
 * Rotates a quaternion by the given angle about the Z axis
 *
 * @param {quat} out quat receiving operation result
 * @param {quat} a quat to rotate
 * @param {number} rad angle (in radians) to rotate
 * @returns {quat} out
 */
quat_p quat_rotate_z(quat_p out, quat_cp a, scalar_t rad) {
	rad *= 0.5;

	scalar_t ax = a[0], ay = a[1], az = a[2], aw = a[3];
	scalar_t bz = sin(rad), bw = cos(rad);

	out[0] = ax * bw + ay * bz;
	out[1] = ay * bw - ax * bz;
	out[2] = az * bw + aw * bz;
	out[3] = aw * bw - az * bz;
	return out;
}

/**
 * Calculates the W component of a quat from the X, Y, and Z components.
 * Assumes that quaternion is 1 unit in length.
 * Any existing W component will be ignored.
 *
 * @param {quat} out the receiving quaternion
 * @param {quat} a quat to calculate W component of
 * @returns {quat} out
 */
quat_p quat_calculate_w(quat_p out, quat_cp a) {
	scalar_t x = a[0], y = a[1], z = a[2];

	out[0] = x;
	out[1] = y;
	out[2] = z;
	out[3] = sqrt(fabs(1.0 - x * x - y * y - z * z));
	return out;
}

/**
 * Performs a spherical linear interpolation between two quat
 *
 * @param {quat} out the receiving quaternion
 * @param {quat} a the first operand
 * @param {quat} b the second operand
 * @param {Number} t interpolation amount between the two inputs
 * @returns {quat} out
 */
quat_p quat_slerp(quat_p out, quat_cp a, quat_cp b, scalar_t t) {
	// benchmarks:
	//    http://jsperf.com/quaternion-slerp-implementations
	scalar_t ax = a[0], ay = a[1], az = a[2], aw = a[3];
	scalar_t bx = b[0], by = b[1], bz = b[2], bw = b[3];

	scalar_t omega, cosom, sinom, scale0, scale1;

	// calc cosine
	cosom = ax * bx + ay * by + az * bz + aw * bw;
	// adjust signs (if necessary)
	if ( cosom < 0.0 ) {
		cosom = -cosom;
		bx = - bx;
		by = - by;
		bz = - bz;
		bw = - bw;
	}
	// calculate coefficients
	if ( (1.0 - cosom) > 0.000001 ) {
		// standard case (slerp)
		omega  = acos(cosom);
		sinom  = sin(omega);
		scale0 = sin((1.0 - t) * omega) / sinom;
		scale1 = sin(t * omega) / sinom;
	} else {
		// "from" and "to" quaternions are very close
		//  ... so we can do a linear interpolation
		scale0 = 1.0 - t;
		scale1 = t;
	}
	// calculate final values
	out[0] = scale0 * ax + scale1 * bx;
	out[1] = scale0 * ay + scale1 * by;
	out[2] = scale0 * az + scale1 * bz;
	out[3] = scale0 * aw + scale1 * bw;

	return out;
}

/**
 * Calculates the inverse of a quat
 *
 * @param {quat} out the receiving quaternion
 * @param {quat} a quat to calculate inverse of
 * @returns {quat} out
 */
quat_p quat_invert(quat_p out, quat_cp a) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
	scalar_t dot = a0*a0 + a1*a1 + a2*a2 + a3*a3;
	scalar_t invDot = dot ? 1.0/dot : 0;

	// TODO: Would be faster to return [0,0,0,0] immediately if dot == 0

	out[0] = -a0*invDot;
	out[1] = -a1*invDot;
	out[2] = -a2*invDot;
	out[3] = a3*invDot;
	return out;
}

/**
 * Calculates the conjugate of a quat
 * If the quaternion is normalized, this function is faster than quat.inverse and produces the same result.
 *
 * @param {quat} out the receiving quaternion
 * @param {quat} a quat to calculate conjugate of
 * @returns {quat} out
 */
quat_p quat_conjugate(quat_p out, quat_cp a) {
	out[0] = -a[0];
	out[1] = -a[1];
	out[2] = -a[2];
	out[3] = a[3];
	return out;
}

/**
 * Creates a quaternion from the given 3x3 rotation matrix.
 *
 * NOTE: The resultant quaternion is not normalized, so you should be sure
 * to renormalize the quaternion yourself where necessary.
 *
 * @param {quat} out the receiving quaternion
 * @param {mat3} m rotation matrix
 * @returns {quat} out
 * @function
 */
quat_p quat_from_mat3(quat_p out, mat3_cp m) {
	// Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
	// article "Quaternion Calculus and Fast Animation".
	scalar_t fTrace = m[0] + m[4] + m[8];
	scalar_t fRoot;

	if ( fTrace > 0.0 ) {
		// |w| > 1/2, may as well choose w > 1/2
		fRoot = sqrt(fTrace + 1.0);  // 2w
		out[3] = 0.5 * fRoot;
		fRoot = 0.5/fRoot;  // 1/(4w)
		out[0] = (m[5]-m[7])*fRoot;
		out[1] = (m[6]-m[2])*fRoot;
		out[2] = (m[1]-m[3])*fRoot;
	} else {
		// |w| <= 1/2
		int i = 0;
		if ( m[4] > m[0] )
			i = 1;
		if ( m[8] > m[i*3+i] )
			i = 2;
		int j = (i+1)%3;
		int k = (i+2)%3;

		fRoot = sqrt(m[i*3+i]-m[j*3+j]-m[k*3+k] + 1.0);
		out[i] = 0.5 * fRoot;
		fRoot = 0.5 / fRoot;
		out[3] = (m[j*3+k] - m[k*3+j]) * fRoot;
		out[j] = (m[j*3+i] + m[i*3+j]) * fRoot;
		out[k] = (m[k*3+i] + m[i*3+k]) * fRoot;
	}

	return out;
}

/**
 * Creates a quaternion from the given euler angle x, y, z.
 *
 * @param {quat} out the receiving quaternion
 * @param {x} Angle to rotate around X axis in degrees.
 * @param {y} Angle to rotate around Y axis in degrees.
 * @param {z} Angle to rotate around Z axis in degrees.
 * @returns {quat} out
 * @function
 */
quat_p quat_from_euler(quat_p out, scalar_t x, scalar_t y, scalar_t z) {
	scalar_t halfToRad = 0.5 * GLM_PI / 180.0;
	x *= halfToRad;
	y *= halfToRad;
	z *= halfToRad;

	scalar_t sx = sin(x);
	scalar_t cx = cos(x);
	scalar_t sy = sin(y);
	scalar_t cy = cos(y);
	scalar_t sz = sin(z);
	scalar_t cz = cos(z);

	out[0] = sx * cy * cz - cx * sy * sz;
	out[1] = cx * sy * cz + sx * cy * sz;
	out[2] = cx * cy * sz - sx * sy * cz;
	out[3] = cx * cy * cz + sx * sy * sz;

	return out;
}

/**
 * Sets a quaternion to represent the shortest rotation from one
 * vector to another.
 *
 * Both vectors are assumed to be unit length.
 *
 * @param {quat} out the receiving quaternion.
 * @param {vec3} a the initial vector
 * @param {vec3} b the destination vector
 * @returns {quat} out
 */
quat_p quat_rotation_to(quat_p out, vec3_cp a, vec3_cp b) {
	scalar_t tmpvec3[3];
	scalar_t xUnitVec3[3] = { 1, 0, 0 };
	scalar_t yUnitVec3[3] = { 0, 1, 0 };

	scalar_t dot = vec3_dot(a, b);
	if (dot < -0.999999) {
		vec3_cross(tmpvec3, xUnitVec3, a);
		if (vec3_len(tmpvec3) < 0.000001)
			vec3_cross(tmpvec3, yUnitVec3, a);
		vec3_normalize(tmpvec3, tmpvec3);
		quat_set_axis_angle(out, tmpvec3, GLM_PI);
		return out;
	} else if (dot > 0.999999) {
		out[0] = 0;
		out[1] = 0;
		out[2] = 0;
		out[3] = 1;
		return out;
	} else {
		vec3_cross(tmpvec3, a, b);
		out[0] = tmpvec3[0];
		out[1] = tmpvec3[1];
		out[2] = tmpvec3[2];
		out[3] = 1 + dot;
		return vec4_normalize(out, out);
	}
}

/**
 * Performs a spherical linear interpolation with two control points
 *
 * @param {quat} out the receiving quaternion
 * @param {quat} a the first operand
 * @param {quat} b the second operand
 * @param {quat} c the third operand
 * @param {quat} d the fourth operand
 * @param {Number} t interpolation amount
 * @returns {quat} out
 */
quat_p quat_sqlerp(quat_p out, quat_cp a, quat_cp b, quat_cp c, quat_cp d, scalar_t t) {
	scalar_t temp1[4];
	scalar_t temp2[4];
	quat_slerp(temp1, a, d, t);
	quat_slerp(temp2, b, c, t);
	quat_slerp(out, temp1, temp2, 2 * t * (1 - t));

	return out;
}

/**
 * Sets the specified quaternion with values corresponding to the given
 * axes. Each axis is a vec3 and is expected to be unit length and
 * perpendicular to all other specified axes.
 *
 * @param {vec3} view  the vector representing the viewing direction
 * @param {vec3} right the vector representing the local "right" direction
 * @param {vec3} up    the vector representing the local "up" direction
 * @returns {quat} out
 */
quat_p quat_set_axes(quat_p out, vec3_cp view, vec3_cp right, vec3_cp up) {
	scalar_t matr[9];
	matr[0] = right[0];
	matr[3] = right[1];
	matr[6] = right[2];

	matr[1] = up[0];
	matr[4] = up[1];
	matr[7] = up[2];

	matr[2] = -view[0];
	matr[5] = -view[1];
	matr[8] = -view[2];

	return vec4_normalize(out, quat_from_mat3(out, matr));
}
