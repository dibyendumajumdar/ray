/* gbtOpticalConstants.h --- constants used in more than one package
   D.Wells, NRAO-CV
   1996-11-14: initial version
   1997-04-29: changed name of include
   1997-07-10: change DTR from constant to expression
   1999-02-03: change rigging angle from 44.0 to 50.8_deg, add BIRDBATH_ANGLE
*/
/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1999 Associated Universities, Inc. Washington DC, USA.

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

#ifndef GBT_OPTICAL_CONSTANTS_H
#define GBT_OPTICAL_CONSTANTS_H

/* TEX{Mathematical and units conversion constants:} */
#define PI 3.1415926535897932
#define DTR (PI / 180.)
#define M2I 39.37
#define I2M (1. / M2I)
#define E3 1000.0
#define E6 1e6

/* TEX{Structural model constants (angles are degrees). The
   \VERB{RIGGING_ANGLE} was set to $50.8^\circ$ or $44^\circ$ here in
   this \VERB{include} file, but is now set to $50.0^\circ$ inside the
   FEM files.} */
#define BIRDBATH_ANGLE 66.0

/* Thermal expansion/contraction constants (RIGGING_TEMP=70_degF) */
/* (LKing 1998-09-24: steel 6.5e-6 in/in/degF, aluminum 12.5E-6 in/in/degF) */
#define RIGGING_TEMP 21.1
/* #define STEEL         1.2e-5  TEMPORARY KLUGE FOR DEBUGGING: */
#define STEEL 0.
#define STAINLESS 1.4e-5
/* #define ALUMINUM      2.3e-5 NOTE SPECIAL KLUGE DEFINITION FOR DEBUGGING! */
#define ALUMINUM STEEL

/* GBT optics parameters (meters & degrees): */
#define GBT_PRIME_FL 60.0
#define GBT_GREG_FL 11.0
#define GBT_EPS 0.528
#define GBT_BETA_ANGLE 5.570
#define GBT_ALPHA_ANGLE 17.899
#define GBT_PF_AXIS_TLT 45.7
#define GBT_SR_AXIS_TLT 36.7

/* TEX{COMSAT engineer J.~Gurney aligned the Gregorian optics of the
   GBT \cite{g00j}; he chose to offset the ``home'' position of the
   subreflector in order to center the range of travel of the
   actuators. The offsets [inches} are:} */
#define X_GURNEY +1.91
#define Y_GURNEY -2.59
#define X_OFFSET -0.95
#define Y_OFFSET -0.50

#endif /* GBT_OPTICAL_CONSTANTS_H */
