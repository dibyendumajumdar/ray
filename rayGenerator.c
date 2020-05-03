/* rayGenerator() --- create lists of rays for plane or spherical waves       
   D.Wells, NRAO-CV
   1993-11: initial version
   1995-07: mods
   1996-03-05: more mods
   1997-01-17: MATRIX_EULER --> VNR_2_MATRIX
   1997-11-21: VNR_2_MATRIX now macro which invokes ANGLES_2_MATRIX
   1998-06-11: 'C3'-->'VC', 'V3'-->'VS', 'MV'-->'MS'
   1998-10-23: remove redundant 'l' in 'g' and 'f' formats
   */

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1996 Associated Universities, Inc. Washington DC, USA.

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

#include "ray.h"
     
/* The following definition is for backward compatibiiity: */
/* Compute rotation matrix AA[3][3] from Euler(?) angles. See "VNR Concise
   Encylopedia of Mathematics'., p.535: E[0] = phi = first Euler angle
   which rotates about the X-axis, E[1]= psi = second Euler angle,
   about Y-axis, E[2] = chi = third angle about Z-axis. Positive
   angles rotate in right- hand-rule sense: if thumb points to +X,
   fingers curl toward positive phi. NOTE: the Euler(?) angle convention
   used here, as given by VNR, differs from that given by Goldstein
   and used in EULER_2_MATRIX(); eventually this difference of
   conventions must be explained properly. */
#define VNR_2_MATRIX(AA, E) \
        ANGLES_2_MATRIX(MC, AA, X,VC(E,0), Y,VC(E,1), Z,VC(E,2), \
                                                        NOP,NOP, NOP,NOP);

#define COLOR_INIT -99999

static int ColorCode = COLOR_INIT;

void rayGenerator(                           /* no value returned            */
		  struct Node *RayBundleSet, /* appends this list of Bundles */
		  char bundle_name[],        /* default used if NULL         */
		  char wave_type[],          /* "plane" || "spherical"       */
		  double wave_point[],       /*         XYZ                  */
		  double wave_direction[],   /*      unit_vector             */
		  double wave_radius,        /* linear  || angular           */
		  double case_step,          /* angular || linear            */
		  int    case_steps,         /*    #steps off-axis           */
		  int    axis_mask,          /* 7=XYZ,4=X,2=Y,6=XY,..        */
		  int    ray_steps,          /*    #rays off-axis            */
		  double taper_angle,        /* of feed horn (radians)       */
		  double taper_db,           /* down at taper_angle          */
		  int ColorCode_1,           /*    first ColorCode           */
		  int ColorCode_2,           /*    last  ColorCode           */
		  enum ColorType assign_by)  /* COLOR_BUNDLE|COLOR_RAY       */
{
    const double tiny=1e-8;
    int i, j, k, l, ir, jr, si, sj, sk; 
    double dl[3], xyz[3], XYZ[3], rad, sum, q, rb, a[3][3], euler[3], edge;
    char temp[NAMEMAX];
    struct Ray *r;
    struct Node *RayBundle;
    MVM_PRIVATE_VARIABLES;
    
    if (RayBundleSet == NULL) {
	printf("rayGenerator: pointer to list of ray bundles is NULL!?!\n");
	exit(EXIT_FAILURE);
    }
    listNodeCheck(RayBundleSet, "r.G: BundleSet");
    if (strlen(RayBundleSet->item) == 0) {
	listFree((char *)RayBundleSet->item, "r.G: empty bundleset string");
	sprintf(temp, "%s, r=%g", wave_type, wave_radius);
	RayBundleSet->item = (char *)listAlloc(strlen(temp)+1, sizeof(char));
	strcpy(RayBundleSet->item, temp);
    }

    if (ColorCode == COLOR_INIT) ColorCode = ColorCode_1;
    q = taper_db / (20.0 * log10(cos(taper_angle)));
    edge = pow(10,(taper_db / 20.0));
    /* if (DEBUG)
       printf("# taper_db=%gdb(%g===>%g), taper_angle=%grad, q = %g\n",
       taper_db, edge, pow(edge, 2.0), taper_angle, q); */
    
    if (strcmp(wave_type, "spherical") == 0) { /* make diverging cones of rays: */
	si = ((axis_mask & 4) != 0);
	for (i = -case_steps*si; i <= +case_steps*si; i++) {
	    sj = ((axis_mask & 2) != 0);
	    for (j = -case_steps*sj; j <= +case_steps*sj; j++) {
		sk = ((axis_mask & 1) != 0);
		for (k = -case_steps*sk; k <= +case_steps*sk; k++) {
		    xyz[0] = (i * case_step);
		    xyz[1] = (j * case_step);
		    xyz[2] = (k * case_step);
		    for (l = 0, rad = 0.0; l < 3; l++) rad += (xyz[l]*xyz[l]);
		    rad = sqrt(rad);
		    if (rad <= ((case_step * case_steps) + tiny)) {
			if (bundle_name != NULL)
			    strcpy(temp, bundle_name);
			else
			    sprintf(temp, "%.4g %.4g %.4g",
				    xyz[0], xyz[1], xyz[2]);
			RayBundle = listInitialize(temp);
			listAppend(RayBundle, RayBundleSet);
			for (l = 0; l < 3; l++) xyz[l] += wave_point[l];
			for (ir = -ray_steps; ir <= +ray_steps; ir++) {
			    dl[0] = 0.0;
			    dl[1] = ir * (wave_radius
					  / ((ray_steps == 0) ? 1 : ray_steps));
			    for (jr = -ray_steps; jr <= +ray_steps; jr++) {
				dl[2] = jr * (wave_radius
					  / ((ray_steps == 0) ? 1 : ray_steps));
				if ((rb = sqrt(dl[1]*dl[1] + dl[2]*dl[2]))
				    > wave_radius) continue;
				r = (struct Ray *)listAlloc(1,sizeof(struct Ray));
				for (l = 0; l < 3; l++) r->T[l] = xyz[l];
				XYZ[0] = 1.0;          /* unit vector */
				XYZ[1] = XYZ[2] = 0.0;
				euler[0] = euler[1] = 0.0; euler[2] = rb;
				VNR_2_MATRIX(a, euler);
				                /* ray tilted off-axis: */
				MATRIX_VECTOR_MULT(VC, r->Q, MC, a, VC, XYZ);  
				euler[0] = (dl[1] == 0.0 &&
					    dl[2] == 0.0) ? 0.0 :
						atan2(dl[2], dl[1]);
				euler[1] = euler[2] = 0.0;
				VNR_2_MATRIX(a, euler); 
				                /* rotated about cone: */
				MATRIX_VECTOR_MULT(VC, XYZ, MC, a, VC, r->Q);
				euler[0] = 0.0;
				euler[1] = -atan2(wave_direction[2],
						  wave_direction[0]);
				euler[2] = -atan2(wave_direction[1],
						  wave_direction[0]);
				VNR_2_MATRIX(a, euler); 
				                /* cone->wave_direction: */
				MATRIX_VECTOR_MULT(VC, r->Q, MC, a, VC, XYZ);
				r->Length = 0.0;
				r->Intensity =
				    pow(cos((rb / wave_radius) * taper_angle),
					2.0 * q);
				r->ColorCode = ColorCode;
				if (assign_by == COLOR_RAY)
				    if (++ColorCode > ColorCode_2)
					ColorCode = ColorCode_1;
				listAppend(r, RayBundle);
			    }
			}
			if (assign_by == COLOR_BUNDLE)
			    if (++ColorCode > ColorCode_2)
				ColorCode = ColorCode_1;
		    }
		}
	    }
	}
    } else if (strcmp(wave_type, "plane") == 0) { /* make plane waves: */
	si = ((axis_mask & 2) != 0);
	for (i = -case_steps*si; i <= +case_steps*si; i++) {
	    XYZ[1] = sin(i * case_step);
	    sj = ((axis_mask & 1) != 0);
	    for (j = -case_steps*sj; j <= +case_steps*sj; j++) {
		XYZ[2] = sin(j * case_step);
		if (sqrt(XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2])
		    < ((case_step * case_steps) + tiny)) {
		    for (l = 0, sum = 0.0, XYZ[0] = 0.0; l < 3; l++) {
			XYZ[l] += wave_direction[l];
			sum += XYZ[l]*XYZ[l];
		    }
		    sum = sqrt(sum);
		    for (l = 0; l < 3; l++) XYZ[l] /= sum;
		    if (bundle_name != NULL)
			strcpy(temp, bundle_name);
		    else
			sprintf(temp, "%8.5f %8.5f", XYZ[1], XYZ[2]);
		    RayBundle = listInitialize(temp);
		    listAppend(RayBundle, RayBundleSet);
		    dl[0] = 0.0;
		    for (ir = -ray_steps; ir <= +ray_steps; ir++) {
			dl[1] = ir * (wave_radius /
				      (ray_steps ? ray_steps : 1));
			for (jr = -ray_steps; jr <= +ray_steps; jr++) {
			    dl[2] = jr * (wave_radius /
					  (ray_steps ? ray_steps : 1));
			    if ((rb = sqrt(dl[1]*dl[1] + dl[2]*dl[2]))
				> (wave_radius + tiny)) continue;
			    r = (struct Ray *) listAlloc(1, sizeof(struct Ray));
			    for (l = 0; l < 3; l++) {
				r->T[l] = wave_point[l] + dl[l];
				r->Q[l] = XYZ[l];
			    }
			    r->Length = 0.0;
			    r->Intensity =
				pow(cos((rb / wave_radius) * taper_angle),
				    2.0 * q);
			    r->ColorCode = ColorCode;
			    if (assign_by == COLOR_RAY)
				if (++ColorCode > ColorCode_2)
				    ColorCode = ColorCode_1;
			    listAppend(r, RayBundle);
			}
		    }
		    if (assign_by == COLOR_BUNDLE)
			if (++ColorCode > ColorCode_2)
			    ColorCode = ColorCode_1;
		}
	    }
	}
    } else {
	printf("rayGenerator: unrecognized wave_type \"%s\"\n",
	       wave_type);
	exit(EXIT_FAILURE);
    }
}
