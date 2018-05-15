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
 * 2x2 Matrix
 * @module mat2
 */

/**
 * Copy the values from one mat2 to another
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the source matrix
 * @returns {mat2} out
 */
mat2_p mat2_copy(mat2_p out, mat2_cp a) {
	out[0] = a[0];
	out[1] = a[1];
	out[2] = a[2];
	out[3] = a[3];
	return out;
}

/**
 * Set a mat2 to the identity matrix
 *
 * @param {mat2} out the receiving matrix
 * @returns {mat2} out
 */
mat2_p mat2_identity(mat2_p out) {
	out[0] = 1;
	out[1] = 0;
	out[2] = 0;
	out[3] = 1;
	return out;
}

/**
 * Set the components of a mat2 to the given values
 *
 * @param {mat2} out the receiving matrix
 * @param {Number} m00 Component in column 0, row 0 position (index 0)
 * @param {Number} m01 Component in column 0, row 1 position (index 1)
 * @param {Number} m10 Component in column 1, row 0 position (index 2)
 * @param {Number} m11 Component in column 1, row 1 position (index 3)
 * @returns {mat2} out
 */
mat2_p mat2_set(mat2_p out, scalar_t m00, scalar_t m01, scalar_t m10, scalar_t m11) {
	out[0] = m00;
	out[1] = m01;
	out[2] = m10;
	out[3] = m11;
	return out;
}

/**
 * Transpose the values of a mat2
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the source matrix
 * @returns {mat2} out
 */
mat2_p mat2_transpose(mat2_p out, mat2_cp a) {
	// If we are transposing ourselves we can skip a few steps but have to cache
	// some values
	if (out == a) {
		scalar_t a1 = a[1];
		out[1] = a[2];
		out[2] = a1;
	} else {
		out[0] = a[0];
		out[1] = a[2];
		out[2] = a[1];
		out[3] = a[3];
	}

	return out;
}

/**
 * Inverts a mat2
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the source matrix
 * @returns {mat2} out
 */
mat2_p mat2_invert(mat2_p out, mat2_cp a) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];

	// Calculate the determinant
	scalar_t det = a0 * a3 - a2 * a1;

	if (!det) {
		return (mat2_p)0;
	}
	det = 1.0 / det;

	out[0] =  a3 * det;
	out[1] = -a1 * det;
	out[2] = -a2 * det;
	out[3] =  a0 * det;

	return out;
}

/**
 * Calculates the adjugate of a mat2
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the source matrix
 * @returns {mat2} out
 */
mat2_p mat2_adjoint(mat2_p out, mat2_cp a) {
	// Caching this value is nessecary if out == a
	scalar_t a0 = a[0];
	out[0] =  a[3];
	out[1] = -a[1];
	out[2] = -a[2];
	out[3] =  a0;

	return out;
}

/**
 * Calculates the determinant of a mat2
 *
 * @param {mat2} a the source matrix
 * @returns {Number} determinant of a
 */
scalar_t mat2_determinant(mat2_cp a) {
	return a[0] * a[3] - a[2] * a[1];
}

/**
 * Multiplies two mat2's
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the first operand
 * @param {mat2} b the second operand
 * @returns {mat2} out
 */
mat2_p mat2_multiply(mat2_p out, mat2_cp a, mat2_cp b) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
	scalar_t b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
	out[0] = a0 * b0 + a2 * b1;
	out[1] = a1 * b0 + a3 * b1;
	out[2] = a0 * b2 + a2 * b3;
	out[3] = a1 * b2 + a3 * b3;
	return out;
}

/**
 * Rotates a mat2 by the given angle
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the matrix to rotate
 * @param {Number} rad the angle to rotate the matrix by
 * @returns {mat2} out
 */
mat2_p mat2_rotate(mat2_p out, mat2_cp a, scalar_t rad) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
	scalar_t s = sin(rad);
	scalar_t c = cos(rad);
	out[0] = a0 *  c + a2 * s;
	out[1] = a1 *  c + a3 * s;
	out[2] = a0 * -s + a2 * c;
	out[3] = a1 * -s + a3 * c;
	return out;
}

/**
 * Scales the mat2 by the dimensions in the given vec2
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the matrix to rotate
 * @param {vec2} v the vec2 to scale the matrix by
 * @returns {mat2} out
 **/
mat2_p mat2_scale(mat2_p out, mat2_cp a, vec2_cp v) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
	scalar_t v0 = v[0], v1 = v[1];
	out[0] = a0 * v0;
	out[1] = a1 * v0;
	out[2] = a2 * v1;
	out[3] = a3 * v1;
	return out;
}

/**
 * Creates a matrix from a given angle
 * This is equivalent to (but much faster than):
 *
 *     mat2.identity(dest);
 *     mat2.rotate(dest, dest, rad);
 *
 * @param {mat2} out mat2 receiving operation result
 * @param {Number} rad the angle to rotate the matrix by
 * @returns {mat2} out
 */
mat2_p mat2_rotation(mat2_p out, scalar_t rad) {
	scalar_t s = sin(rad);
	scalar_t c = cos(rad);
	out[0] = c;
	out[1] = s;
	out[2] = -s;
	out[3] = c;
	return out;
}

/**
 * Creates a matrix from a vector scaling
 * This is equivalent to (but much faster than):
 *
 *     mat2.identity(dest);
 *     mat2.scale(dest, dest, vec);
 *
 * @param {mat2} out mat2 receiving operation result
 * @param {vec2} v Scaling vector
 * @returns {mat2} out
 */
mat2_p mat2_scaling(mat2_p out, vec2_cp v) {
	out[0] = v[0];
	out[1] = 0;
	out[2] = 0;
	out[3] = v[1];
	return out;
}

/**
 * Returns Frobenius norm of a mat2
 *
 * @param {mat2} a the matrix to calculate Frobenius norm of
 * @returns {Number} Frobenius norm
 */
scalar_t mat2_frob(mat2_cp a) {
	return(sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2) + pow(a[3], 2)));
}

/**
 * Returns L, D and U matrices (Lower triangular, Diagonal and Upper triangular) by factorizing the input matrix
 * @param {mat2} L the lower triangular matrix
 * @param {mat2} D the diagonal matrix
 * @param {mat2} U the upper triangular matrix
 * @param {mat2} a the input matrix to factorize
 */

void mat2_LDU(mat2_p L, mat2_p D, mat2_p U, mat2_cp a) {
	L[2] = a[2]/a[0];
	U[0] = a[0];
	U[1] = a[1];
	U[3] = a[3] - L[2] * U[1];
}

/**
 * Adds two mat2's
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the first operand
 * @param {mat2} b the second operand
 * @returns {mat2} out
 */
mat2_p mat2_add(mat2_p out, mat2_cp a, mat2_cp b) {
	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
	out[2] = a[2] + b[2];
	out[3] = a[3] + b[3];
	return out;
}

/**
 * Subtracts matrix b from matrix a
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the first operand
 * @param {mat2} b the second operand
 * @returns {mat2} out
 */
mat2_p mat2_sub(mat2_p out, mat2_cp a, mat2_cp b) {
	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
	out[2] = a[2] - b[2];
	out[3] = a[3] - b[3];
	return out;
}

/**
 * Returns whether or not the matrices have exactly the same elements in the same position (when compared with ===)
 *
 * @param {mat2} a The first matrix.
 * @param {mat2} b The second matrix.
 * @returns {Boolean} True if the matrices are equal, false otherwise.
 */
int mat2_exact_equals(mat2_cp a, mat2_cp b) {
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

/**
 * Returns whether or not the matrices have approximately the same elements in the same position.
 *
 * @param {mat2} a The first matrix.
 * @param {mat2} b The second matrix.
 * @returns {Boolean} True if the matrices are equal, false otherwise.
 */
int mat2_equals(mat2_cp a, mat2_cp b) {
	scalar_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
	scalar_t b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
	return (fabs(a0 - b0) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a0), fabs(b0))) &&
		fabs(a1 - b1) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a1), fabs(b1))) &&
		fabs(a2 - b2) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a2), fabs(b2))) &&
		fabs(a3 - b3) <= GLM_EPSILON*fmax(1.0, fmax(fabs(a3), fabs(b3))));
}

/**
 * Multiply each element of the matrix by a scalar.
 *
 * @param {mat2} out the receiving matrix
 * @param {mat2} a the matrix to scale
 * @param {Number} b amount to scale the matrix's elements by
 * @returns {mat2} out
 */
mat2_p mat2_multiply_scalar(mat2_p out, mat2_cp a, scalar_t b) {
	out[0] = a[0] * b;
	out[1] = a[1] * b;
	out[2] = a[2] * b;
	out[3] = a[3] * b;
	return out;
}

/**
 * Adds two mat2's after multiplying each element of the second operand by a scalar value.
 *
 * @param {mat2} out the receiving vector
 * @param {mat2} a the first operand
 * @param {mat2} b the second operand
 * @param {Number} scale the amount to scale b's elements by before adding
 * @returns {mat2} out
 */
mat2_p mat2_madd(mat2_p out, mat2_cp a, mat2_cp b, scalar_t scale) {
	out[0] = a[0] + (b[0] * scale);
	out[1] = a[1] + (b[1] * scale);
	out[2] = a[2] + (b[2] * scale);
	out[3] = a[3] + (b[3] * scale);
	return out;
}
