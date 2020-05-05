/* mathVectorMatrix.h --- C macros for 3-vector and 3x3-matrix operations

   1996-11-14: D.Wells, NRAO-CV, early version
   1997-01-10: modify small-angle algorithm
   1997-01-17: re-install old VNR Euler-angle algorithm as kluge
   1997-02-06: added MATRIX_DETERMINANT()
   1997-05-03: many name changes
   1997-06-13: added <stdio.h> for printf() use, "%12.6lg"--->"%12.6g"
   1997-07-11: change indexed arrays to "subscripted" scalars
   1997-07-23: added MATRIX_MULT_9, implemented C3 & V3 vector "types"
   1997-07-24: MATRIX_ROTATE-->MATRIX_VECTOR_MULT, added MATRIX_TRANSPOSE, etc
   1997-10-27: added MC & MV matrix "types", change MATRIX_COPY
   1997-10-28: changed matrix init and multiply code to new macros,
	       timing actually slowed slightly. sigh...
   1997-11-21: added matrix type codes to more matrix macros
   1998-05-15: final semicolon removed from VECTOR3() and MATRIX3(),
	       MATRIX_2_ARRAY() and EULER_2_MATRIX_OLD() deleted
   1998-05-29: 'C3'-->'VC', 'V3'-->'VS', 'MV'-->'MS',
	       'VECTOR3'-->'VECTOR_VS', 'MATRIX3'-->'MATRIX_MS'
   1998-06-10: MATRIX_MULT() uses temp matrix so arg C can be same as A|B.
   1999-08-07: VECTOR_ANGLE scale ratio (1-1e-12) so collinear not yield NaN
   2001-03-21: moved vector/matrix debug macros to gbtDebugMacros.h
   2003-01-31: delete ## from macros in places where gcc3.2 gives warnings
*/

#ifndef VECTOR_MATRIX_H
#define VECTOR_MATRIX_H

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 2003 Associated Universities, Inc. Washington DC, USA.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   Correspondence concerning GBT software should be addressed as follows:
   GBT Operations, National Radio Astronomy Observatory, P. O. Box 2,
   Green Bank, WV 24944-0002 USA

   -=-=-=-=-=-=-=-=-=- End Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=- */

#include "gbtOpticalConstants.h"
#include <math.h>
#include <stdio.h>

#define VECTOR_VS(V) double V##_0, V##_1, V##_2

#define MATRIX_MS(M) double M##_00, M##_01, M##_02, M##_10, M##_11, M##_12, M##_20, M##_21, M##_22

#define MVM_PRIVATE_VARIABLES                                                                                          \
	struct {                                                                                                       \
		double vl;                                                                                             \
		VECTOR_VS(vtemp);                                                                                      \
		MATRIX_MS(mtemp);                                                                                      \
		double xca, xsa;                                                                                       \
		MATRIX_MS(xt2);                                                                                        \
		double yca, ysa;                                                                                       \
		MATRIX_MS(yt2);                                                                                        \
		double zca, zsa;                                                                                       \
		MATRIX_MS(zt2);                                                                                        \
	} MVM

/* define the two types of vector subscripts: */
#define VC(V, S) V[S]
#define VS(V, S) V##_##S

/* define the two types of matrix subscripts: */
#define MC(M, I, J) M[I][J]
#define MS(M, I, J) M##_##I##J

#define VECTOR_INIT(TV, V, v0, v1, v2)                                                                                 \
	TV(V, 0) = v0;                                                                                                 \
	TV(V, 1) = v1;                                                                                                 \
	TV(V, 2) = v2;

#define VECTOR_FILL(TV, V, v0) TV(V, 0) = TV(V, 1) = TV(V, 2) = v0;

#define VECTOR_ADD(TV, V, TA, A, TB, B)                                                                                \
	TV(V, 0) = TA(A, 0) + TB(B, 0);                                                                                \
	TV(V, 1) = TA(A, 1) + TB(B, 1);                                                                                \
	TV(V, 2) = TA(A, 2) + TB(B, 2);

#define VECTOR_SUB(TV, V, TA, A, TB, B)                                                                                \
	TV(V, 0) = TA(A, 0) - TB(B, 0);                                                                                \
	TV(V, 1) = TA(A, 1) - TB(B, 1);                                                                                \
	TV(V, 2) = TA(A, 2) - TB(B, 2);

#define VECTOR_DOT(TA, A, TB, B) (TA(A, 0) * TB(B, 0) + TA(A, 1) * TB(B, 1) + TA(A, 2) * TB(B, 2))

#define VECTOR_SCALE(TV, V, s, TA, A)                                                                                  \
	TV(V, 0) = s * TA(A, 0);                                                                                       \
	TV(V, 1) = s * TA(A, 1);                                                                                       \
	TV(V, 2) = s * TA(A, 2);

/* add vector A[3] times a constant to vector V[3] */
#define VECTOR_EXTEND(TV, V, TA, A, c)                                                                                 \
	TV(V, 0) += c * TA(A, 0);                                                                                      \
	TV(V, 1) += c * TA(A, 1);                                                                                      \
	TV(V, 2) += c * TA(A, 2);

#define VECTOR_LENGTH(TV, V) (sqrt(VECTOR_DOT(TV, V, TV, V)))

#define VECTOR_UNIT(TA, A, TB, B)                                                                                      \
	MVM.vl = (1. / VECTOR_LENGTH(TB, B));                                                                          \
	VECTOR_SCALE(TA, A, MVM.vl, TB, B)

#define VECTOR_ANGLE(TA, A, TB, B)                                                                                     \
	(acos((VECTOR_DOT(TA, A, TB, B) * (1. - 1e-14)) / (VECTOR_LENGTH(TA, A) * VECTOR_LENGTH(TB, B))))

#define VECTOR_CROSS(TV, V, TA, A, TB, B)                                                                              \
	TV(V, 0) = (TA(A, 1) * TB(B, 2) - TA(A, 2) * TB(B, 1));                                                        \
	TV(V, 1) = (TA(A, 2) * TB(B, 0) - TA(A, 0) * TB(B, 2));                                                        \
	TV(V, 2) = (TA(A, 0) * TB(B, 1) - TA(A, 1) * TB(B, 0));

#define VECTOR_COPY(TV, V, TA, A)                                                                                      \
	TV(V, 0) = TA(A, 0);                                                                                           \
	TV(V, 1) = TA(A, 1);                                                                                           \
	TV(V, 2) = TA(A, 2);

#define VECTOR_PERMUTE(TV, V, I, J, K, FI, FJ, FK, TA, A)                                                              \
	VS(MVM.vtemp, 0) = FI * TA(A, I);                                                                              \
	VS(MVM.vtemp, 1) = FJ * TA(A, J);                                                                              \
	VS(MVM.vtemp, 2) = FK * TA(A, K);                                                                              \
	VECTOR_COPY(TV, V, VS, MVM.vtemp);

/* Macro to store 9 numbers into a 3x3 matrix: */
#define MATRIX_INIT(TM, M, M00, M01, M02, M10, M11, M12, M20, M21, M22)                                                \
	TM(M, 0, 0) = M00;                                                                                             \
	TM(M, 0, 1) = M01;                                                                                             \
	TM(M, 0, 2) = M02;                                                                                             \
	TM(M, 1, 0) = M10;                                                                                             \
	TM(M, 1, 1) = M11;                                                                                             \
	TM(M, 1, 2) = M12;                                                                                             \
	TM(M, 2, 0) = M20;                                                                                             \
	TM(M, 2, 1) = M21;                                                                                             \
	TM(M, 2, 2) = M22;

#define MATRIX_COPY(TA, A, TB, B)                                                                                      \
	TA(A, 0, 0) = TB(B, 0, 0);                                                                                     \
	TA(A, 0, 1) = TB(B, 0, 1);                                                                                     \
	TA(A, 0, 2) = TB(B, 0, 2);                                                                                     \
	TA(A, 1, 0) = TB(B, 1, 0);                                                                                     \
	TA(A, 1, 1) = TB(B, 1, 1);                                                                                     \
	TA(A, 1, 2) = TB(B, 1, 2);                                                                                     \
	TA(A, 2, 0) = TB(B, 2, 0);                                                                                     \
	TA(A, 2, 1) = TB(B, 2, 1);                                                                                     \
	TA(A, 2, 2) = TB(B, 2, 2);

#define MATRIX_PERMUTE(TMP, MP, I, J, K, TM, M)                                                                        \
	MS(MVM.mtemp, 0, 0) = TM(M, I, I);                                                                             \
	MS(MVM.mtemp, 0, 1) = TM(M, I, J);                                                                             \
	MS(MVM.mtemp, 0, 2) = TM(M, I, K);                                                                             \
	MS(MVM.mtemp, 1, 0) = TM(M, J, I);                                                                             \
	MS(MVM.mtemp, 1, 1) = TM(M, J, J);                                                                             \
	MS(MVM.mtemp, 1, 2) = TM(M, J, K);                                                                             \
	MS(MVM.mtemp, 2, 0) = TM(M, K, I);                                                                             \
	MS(MVM.mtemp, 2, 1) = TM(M, K, J);                                                                             \
	MS(MVM.mtemp, 2, 2) = TM(M, K, K);                                                                             \
	MATRIX_COPY(TMP, MP, MS, MVM.mtemp);

#define MATRIX_DETERMINANT(TA, A)                                                                                      \
	(TA(A, 0, 0) * TA(A, 1, 1) * TA(A, 2, 2) - TA(A, 0, 0) * TA(A, 1, 2) * TA(A, 2, 1) -                           \
	 TA(A, 0, 1) * TA(A, 1, 0) * TA(A, 2, 2) + TA(A, 0, 1) * TA(A, 1, 2) * TA(A, 2, 0) +                           \
	 TA(A, 0, 2) * TA(A, 1, 0) * TA(A, 2, 1) - TA(A, 0, 2) * TA(A, 1, 1) * TA(A, 2, 0))

#define ROTATION_ANGLE(TRM, RM) (acos(0.5 * (TRM(RM, 0, 0) + TRM(RM, 1, 1) + TRM(RM, 2, 2) - 1.0)))

/* Macros to multiply A[3][3] by B[3][3] and return result in C[][]:
   for (m = 0; m < 3; m++) for (n = 0; n < 3; n++)
       C[n][m] = A[n][0*B[0][m] + A[n][1]*B[1][m] + A[n][2]*B[2][m]; */
#define MATRIX_MULT(TC, C, TA, A, TB, B)                                                                               \
	MS(MVM.mtemp, 0, 0) = TA(A, 0, 0) * TB(B, 0, 0) + TA(A, 0, 1) * TB(B, 1, 0) + TA(A, 0, 2) * TB(B, 2, 0);       \
	MS(MVM.mtemp, 1, 0) = TA(A, 1, 0) * TB(B, 0, 0) + TA(A, 1, 1) * TB(B, 1, 0) + TA(A, 1, 2) * TB(B, 2, 0);       \
	MS(MVM.mtemp, 2, 0) = TA(A, 2, 0) * TB(B, 0, 0) + TA(A, 2, 1) * TB(B, 1, 0) + TA(A, 2, 2) * TB(B, 2, 0);       \
	MS(MVM.mtemp, 0, 1) = TA(A, 0, 0) * TB(B, 0, 1) + TA(A, 0, 1) * TB(B, 1, 1) + TA(A, 0, 2) * TB(B, 2, 1);       \
	MS(MVM.mtemp, 1, 1) = TA(A, 1, 0) * TB(B, 0, 1) + TA(A, 1, 1) * TB(B, 1, 1) + TA(A, 1, 2) * TB(B, 2, 1);       \
	MS(MVM.mtemp, 2, 1) = TA(A, 2, 0) * TB(B, 0, 1) + TA(A, 2, 1) * TB(B, 1, 1) + TA(A, 2, 2) * TB(B, 2, 1);       \
	MS(MVM.mtemp, 0, 2) = TA(A, 0, 0) * TB(B, 0, 2) + TA(A, 0, 1) * TB(B, 1, 2) + TA(A, 0, 2) * TB(B, 2, 2);       \
	MS(MVM.mtemp, 1, 2) = TA(A, 1, 0) * TB(B, 0, 2) + TA(A, 1, 1) * TB(B, 1, 2) + TA(A, 1, 2) * TB(B, 2, 2);       \
	MS(MVM.mtemp, 2, 2) = TA(A, 2, 0) * TB(B, 0, 2) + TA(A, 2, 1) * TB(B, 1, 2) + TA(A, 2, 2) * TB(B, 2, 2);       \
	MATRIX_COPY(TC, C, MS, MVM.mtemp);

#define MATRIX_MULT_9(TC, C, TA, A, m00, m01, m02, m10, m11, m12, m20, m21, m22)                                       \
	TC(C, 0, 0) = TA(A, 0, 0) * m00 + TA(A, 0, 1) * m10 + TA(A, 0, 2) * m20;                                       \
	TC(C, 1, 0) = TA(A, 1, 0) * m00 + TA(A, 1, 1) * m10 + TA(A, 1, 2) * m20;                                       \
	TC(C, 2, 0) = TA(A, 2, 0) * m00 + TA(A, 2, 1) * m10 + TA(A, 2, 2) * m20;                                       \
	TC(C, 0, 1) = TA(A, 0, 0) * m01 + TA(A, 0, 1) * m11 + TA(A, 0, 2) * m21;                                       \
	TC(C, 1, 1) = TA(A, 1, 0) * m01 + TA(A, 1, 1) * m11 + TA(A, 1, 2) * m21;                                       \
	TC(C, 2, 1) = TA(A, 2, 0) * m01 + TA(A, 2, 1) * m11 + TA(A, 2, 2) * m21;                                       \
	TC(C, 0, 2) = TA(A, 0, 0) * m02 + TA(A, 0, 1) * m12 + TA(A, 0, 2) * m22;                                       \
	TC(C, 1, 2) = TA(A, 1, 0) * m02 + TA(A, 1, 1) * m12 + TA(A, 1, 2) * m22;                                       \
	TC(C, 2, 2) = TA(A, 2, 0) * m02 + TA(A, 2, 1) * m12 + TA(A, 2, 2) * m22;

/* Generate arbitrary 5-angle rotation matrix.
       See Euler-angle macro below for typical call. */
#define ABOUT_X(TRM, RM, A)                                                                                            \
	MVM.xca = cos(A);                                                                                              \
	MVM.xsa = sin(A);                                                                                              \
	MATRIX_MULT_9(MS, MVM.xt2, TRM, RM, 1., 0., 0., 0., MVM.xca, -MVM.xsa, 0., +MVM.xsa, MVM.xca)                  \
	MATRIX_COPY(TRM, RM, MS, MVM.xt2);

#define ABOUT_Y(TRM, RM, A)                                                                                            \
	MVM.yca = cos(A);                                                                                              \
	MVM.ysa = sin(A);                                                                                              \
	MATRIX_MULT_9(MS, MVM.yt2, TRM, RM, MVM.yca, 0., -MVM.ysa, 0., 1., 0., +MVM.ysa, 0., MVM.yca)                  \
	MATRIX_COPY(TRM, RM, MS, MVM.yt2);

#define ABOUT_Z(TRM, RM, A)                                                                                            \
	MVM.zca = cos(A);                                                                                              \
	MVM.zsa = sin(A);                                                                                              \
	MATRIX_MULT_9(MS, MVM.zt2, TRM, RM, MVM.zca, -MVM.zsa, 0., +MVM.zsa, MVM.zca, 0., 0., 0., 1.)                  \
	MATRIX_COPY(TRM, RM, MS, MVM.zt2);
#define ABOUT_NOP(TRM, RM, A)
/* Use An,Rn="NOP,NOP" for unused rotations (see EULER_2_MATRIX). */
#define ANGLES_2_MATRIX(TR, R, A1, R1, A2, R2, A3, R3, A4, R4, A5, R5)                                                 \
	MATRIX_INIT(TR, R, 1., 0., 0., 0., 1., 0., 0., 0., 1.);                                                        \
	ABOUT_##A1(TR, R, R1);                                                                                         \
	ABOUT_##A2(TR, R, R2);                                                                                         \
	ABOUT_##A3(TR, R, R3);                                                                                         \
	ABOUT_##A4(TR, R, R4);                                                                                         \
	ABOUT_##A5(TR, R, R5);

/* Compute rotation matrix AA[3][3] from three Euler angles E[3].
   AA = Z(alpha)X(beta)Z(gamma)
   E[0] = alpha = first rotation matrix, about z-axis, 0 <= alpha < 2*pi
   E[1] = beta = second rotation matrix, about x-axis, 0 <= beta  <=  pi
   E[2] = gamma = third rotation matrix, about Z-axis, 0 <= gamma < 2*pi
   Rotations obey right-hand rule (thumb to +axis, finger curl +rotation). */
#define EULER_2_MATRIX(TAA, AA, TE, E)                                                                                 \
	ANGLES_2_MATRIX(TAA, AA, Z, TE(E, 0), X, TE(E, 1), Z, TE(E, 2), NOP, NOP, NOP, NOP);

/* Compute three Euler angles from rotation matrix AA[3][3]. */
#define LT_EPS(A) (fabs(A) < 1.e-7)
#define MATRIX_2_EULER(TE, E, TAA, AA)                                                                                 \
	TE(E, 1) = acos(TAA(AA, 2, 2));                                                                                \
	if (LT_EPS(TE(E, 1))) {                                                                                        \
		TE(E, 0) = 0.; /* alpha is indeterminate for beta=0. */                                                \
		TE(E, 2) = +atan2(TAA(AA, 1, 0), TAA(AA, 0, 0));                                                       \
	} else {                                                                                                       \
		if (LT_EPS(TE(E, 1) - PI)) {                                                                           \
			TE(E, 0) = 0.; /* also indeterminate if beta=PI. */                                            \
			TE(E, 2) = -atan2(TAA(AA, 1, 0), TAA(AA, 0, 0));                                               \
		} else {                                                                                               \
			TE(E, 0) = (LT_EPS(TAA(AA, 0, 2)) && LT_EPS(TAA(AA, 1, 2)))                                    \
				       ? 0.                                                                            \
				       : atan2(TAA(AA, 0, 2), -TAA(AA, 1, 2));                                         \
			if (TE(E, 0) < 0.)                                                                             \
				TE(E, 0) += 2 * PI;                                                                    \
			TE(E, 2) = (LT_EPS(TAA(AA, 2, 0)) && LT_EPS(TAA(AA, 2, 1)))                                    \
				       ? 0.                                                                            \
				       : atan2(TAA(AA, 2, 0), TAA(AA, 2, 1));                                          \
		}                                                                                                      \
	};                                                                                                             \
	if (TE(E, 2) < 0.)                                                                                             \
		TE(E, 2) += 2 * PI;

/* Compute rotation matrix RM[3][3] from small-angle approximation.
   The errors committed in this approximation will be comparable to the
   squares of the angles in radians. I.e., angles of order 10^-3 radians
   will produce rotation matrices accurate to about 10^-6 radian.
   E[0] about X-axis, E[1] about Y, E[2] about Z, right-hand rule.
			    1.,    +E[2, -E[1,
			    -E[2, 1.,    +E[0,
			    +E[1, -E[0, 1.       */
#define SMALL_2_MATRIX(TRM, RM, TE, E)                                                                                 \
	ANGLES_2_MATRIX(TRM, RM, X, TE(E, 0), Y, TE(E, 1), Z, TE(E, 2), NOP, NOP, NOP, NOP)

/* Compute small-angles from rotation matrix RM[3][3].
   Eq.4-91 in Goldstein's "Classical Mechanics" (1950) can be expressed
   (in the notation used here) as
		   | 0     +E[2  -E[1] |
	     E[] = | -E[2] 0     +E[0] |
		   | +E[1] -E[0] 0     |
   This approximation is valid to the extent that sinE~=0, cosE~=1,
   i.e., to the extent that sin(E[i])*sin(E[j])<<0. */
#define MATRIX_2_SMALL(TE, E, TRM, RM)                                                                                 \
	TE(E, 0) = asin(0.5 * (TRM(RM, 1, 2) - TRM(RM, 2, 1)));                                                        \
	TE(E, 1) = asin(0.5 * (TRM(RM, 2, 0) - TRM(RM, 0, 2)));                                                        \
	TE(E, 2) = asin(0.5 * (TRM(RM, 0, 1) - TRM(RM, 1, 0)));

/* Multiply (rotation) matrix R[3][3] by vector V[3] to get vector
   VP[3]. See VNR p.534 -- convention is that coordinates of V[] are
   (x,y,z) and coordinates of VP[] are (x^*,y^*,z^*). */
#define MATRIX_VECTOR_MULT(TVP, VP, TR, R, TV, V)                                                                      \
	TVP(VP, 0) = TV(V, 0) * TR(R, 0, 0) + TV(V, 1) * TR(R, 1, 0) + TV(V, 2) * TR(R, 2, 0);                         \
	TVP(VP, 1) = TV(V, 0) * TR(R, 0, 1) + TV(V, 1) * TR(R, 1, 1) + TV(V, 2) * TR(R, 2, 1);                         \
	TVP(VP, 2) = TV(V, 0) * TR(R, 0, 2) + TV(V, 1) * TR(R, 1, 2) + TV(V, 2) * TR(R, 2, 2);

#define MATRIX_TRANSPOSE(TMT, MT, TM, M)                                                                               \
	TMT(MT, 0, 0) = TM(M, 0, 0);                                                                                   \
	TMT(MT, 0, 1) = TM(M, 1, 0);                                                                                   \
	TMT(MT, 0, 2) = TM(M, 2, 0);                                                                                   \
	TMT(MT, 1, 0) = TM(M, 0, 1);                                                                                   \
	TMT(MT, 1, 1) = TM(M, 1, 1);                                                                                   \
	TMT(MT, 1, 2) = TM(M, 2, 1);                                                                                   \
	TMT(MT, 2, 0) = TM(M, 0, 2);                                                                                   \
	TMT(MT, 2, 1) = TM(M, 1, 2);                                                                                   \
	TMT(MT, 2, 2) = TM(M, 2, 2);

#endif /* VECTOR_MATRIX_H */
