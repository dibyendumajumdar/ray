/* ray.h -- Include file for the 'ray' (ray-tracing) package.
   D.Wells, NRAO-CV
   1997-05-03: many name changes
   1997-05-30: removed changes to rayGenerator()
   1998-10-25: revised list utilities, new functions.
*/

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1997 Associated Universities, Inc. Washington DC, USA.

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

#ifndef RAY_H
#define RAY_H

#include "mathVectorMatrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
#define NAMEMAX 21
#define FORMMAX 100
#define AMNIMAX 9
enum VignetteType { VIGN_CYLINDER, VIGN_CONE };
enum ColorType { COLOR_BUNDLE, COLOR_RAY };
enum PlotType { ORTHOGRAPHIC, ORTHOGRAPHIC_X, PERSPECTIVE };

struct Surface {
	char name[NAMEMAX]; /* user-supplied descriptive string */
	double c_1;	    /* curvature c=1/r; positive if center of curvature is
			       to right of surface. c_1 for following surface */
	double eps;	    /* numerical eccentricity of the conic section. */
	double a_2;	    /* coefficients of radially-symmetric */
	double a_4;	    /*       aspheric deformation series. */
	double mu_1;	    /* $mu_1\equiv N/N_1$, N_1 is index to right
			       of following surface */
	double S[3];	    /* absolute vector to vertex of this surface. */
	double E[3];	    /* tilt of next vertex, 3 Euler Angles (see
			       "General Ray-Tracing Procedure",
			       G.H.Spencer and M.V.R.K.Murty, JOSA 52,
			       pp.672-678, (June 1962)). */
	enum VignetteType vign_type;
	double vign_origin[3];
	double vign_vector[3];
	double vign_radius;
};

struct Ray {
	double T[3];	  /* $T \equiv (x, y, z)$ is the vector from the
			     vertex of the first surface to the point of
			     incidence of the ray on this surface.          */
	double Q[3];	  /* $Q \equiv (X, Y, Z)$ is the unit vector along
			     the ray to the right of the first surface. X,Y,Z
			     are direction cosines; ray path is $T + l * Q$ */
	double direction; /* experimental, +1->+X, -1-->-X                  */
	double Length;	  /* cumulative path length as ray traverses system */
	double Intensity; /* same as power for radio                        */
	int ColorCode;	  /* color code to use when plotting this ray       */
};

struct Segment {
	double T1[3];  /* x,y,z position of starting point of ray */
	double T2[3];  /* x,y,z position of   ending point of ray */
	int ColorCode; /* color code to use when plotting this ray */
};

struct Focus {
	char name[NAMEMAX]; /* name_string for the ray bundle */
	int n;		    /* number of rays in the bundle   */
	double xyz[3];	    /* position of focus              */
	double xyzs[3];	    /* std_err of focus position      */
	double lc;	    /* mean pathlength to focus       */
	double ls;	    /* spherical wave rms wrt focus   */
};

struct Plane {
	char name[NAMEMAX];   /* name_string for the ray bundle   */
	int n;		      /* number of rays in the bundle     */
	int nu;		      /* number of Seidel terms           */
	double Amni[AMNIMAX]; /* Seidel terms for pathlengths     */
	double Asig[AMNIMAX]; /* Sigmas for Amni[AMNIMAX]         */
	int mni[AMNIMAX];     /* indicies of Seidel terms         */
	double sig;	      /* rms of (pathlength-sum(Amni[]))  */
	double yz_radius;     /* Zernike radius in rayGetPlanes() */
};

#ifdef SEIDEL
struct Seidel_item {
	int mni;    /* A_mni index                     */
	const char *name; /* names of Seidel aberrrations    */
};
static struct Seidel_item Seidel_list[] = {{0000, "Zero_Pt"}, {2000, "Defoc"},	     {4000, "Sph_Ab"},
					   {1010, "Tilt"},    {1011, " <sin>"},	     {3010, "Coma"},
					   {3011, " <sin>"},  {2020, "Astigmatism"}, {2021, " <sin>"}};
static const int Seidel_n = sizeof(Seidel_list) / sizeof(struct Seidel_item);
#endif

#define NODECHECK 229032475
#define NODETEXT 7
#define NODE_DEBUG FALSE

struct Node {
	void *item;	     /* pointer to an "item" struct */
	struct Node *next;   /* pointer to next item in list */
	int nodecheck;	     /* private value to verify struct Node *p */
	char text[NODETEXT]; /* descriptive comment about function of this Node */
};

#define NODE_MAX 2000000
struct {
	int N_CALLOC;
	void *CALLOC_LIST[NODE_MAX]; /* private list of calloc() results */
} NODE_PRIVATE;

/* -=-=-=-=-=-=-=-=-=-=- Function Prototypes: -=-=-=-=-=-=-=-=-=-=-=- */
#ifdef __cplusplus
extern "C" {
#endif /* endif cplusplus */

struct Node *rayTrace(struct Node *RayBundleSet, /* the list of input rays  */
		      struct Node *System,	 /* the list of surfaces    */
		      double tolerance,		 /* for aspheric intercepts */
		      struct Node *Segments)	 /* init this list! appends */
    ;

void rayGenerator(			     /* no value returned            */
		  struct Node *RayBundleSet, /* appends this list of Bundles */
		  char bundle_name[],	     /* default used if NULL         */
		  char wave_type[],	     /* "plane" || "spherical"       */
		  double wave_point[],	     /*         XYZ                  */
		  double wave_direction[],   /*      unit_vector             */
		  double wave_radius,	     /* linear  || angular           */
		  double case_step,	     /* angular || linear            */
		  int case_steps,	     /*    #steps off-axis           */
		  int axis_mask,	     /* 7=XYZ,4=X,2=Y,6=XY,..        */
		  int ray_steps,	     /*    #rays off-axis            */
		  double taper_angle,	     /* of feed horn (radians)       */
		  double taper_db,	     /* down at taper_angle          */
		  int ColorCode_1,	     /*    first ColorCode           */
		  int ColorCode_2,	     /*    last  ColorCode           */
		  enum ColorType assign_by)  /* COLOR_BUNDLE|COLOR_RAY       */
    ;

void rayPrtBundles(struct Node *set, /* list of lists of rays      */
		   int d)	     /* digits after decimal point */
    ;
struct Node *rayAddSurface(struct Node *list,		   /* list of surfaces */
			   char surfname[],		   /* descriptive string */
			   double curvature,		   /* 1/r */
			   double epsilon,		   /* eccentricity of conic */
			   double A_2,			   /* deformation terms */
			   double A_4, double index_ratio, /* N/N_1, -1 means mirror */
			   double S[],			   /* XYZ of vertex */
			   double E[],			   /* Euler angles of vertex tilt */
			   enum VignetteType vign_type,	   /* VIGN_CYLINDER | VIGN_CONE */
			   double VO[],			   /* vignette origin XYZ */
			   double VV[],			   /* vignette direction cosines */
			   double VR)			   /* radius, linear | radians   */
    ;
void rayPrtSystem(struct Node *list, /* list of surfaces           */
		  int d)	     /* digits after decimal point */
    ;
void rayPrtSegments(struct Node *list, /* list of ray segments       */
		    int d)	       /* digits after decimal point */
    ;
struct Node *rayGetFoci(			   /* returns list of Focus */
			struct Node *RayBundleSet, /* list of lists of rays */
			double tolerance,	   /* position accuracy     */
			struct Node *Segments)	   /* appends to Segments   */
    ;
void rayPrtFoci(struct Node *foci_list, /* list of Focus structs      */
		int d)			/* digits after decimal point */
    ;
struct Node *rayGetPlanes(			     /* return list of planes */
			  struct Node *RayBundleSet, /* list of lists of rays */
			  double yz_axis[],	     /* Zernike analysis axis */
			  double yz_radius,	     /* Zernike radius        */
			  char print_mode[])	     /* "verbose"|"silent"    */
    ;
void rayPrtPlanes(struct Node *planes_list, /* list of Plane structs      */
		  int d)		    /* digits after decimal point */
    ;
int rayPltPS(			    /* returns non-zero on error     */
	     struct Node *segments, /* list of line segments         */
	     double height_cm,	    /* height of PS in centimeters   */
	     double width_cm,	    /* width  of PS in centimeters   */
	     double width,	    /* width in units of System      */
	     double to_point[],	    /* XYZ of point in center of PS  */
	     char plt_mode[],	    /* "Orthographic"|"Perspective"  */
	     char psname[])	    /* Postscript output file        */
    ;
void rayPltPSDeleteHues(void) /* no args, no return value */
    ;
void rayPltSystem(			 /* returns lines in segments list */
		  struct Node *system,	 /* list of optical elements       */
		  struct Node *segments) /* list of line segments          */
    ;
struct Node *listNodeCreate(		 /* return ptr to calloc(Node) */
			    const char *text) /* descriptive text re Node   */
    ;
struct Node *listInitialize(		/* return ptr to new list      */
			    const char *name) /* descriptive string for list */
    ;
struct Node *listAppend(void *newitem,	   /* ptr to new item struct */
			struct Node *list) /* ptr to existing list   */
    ;
struct Node *listInsertAfter(void *newitem,  /* ptr to new item struct */
			     struct Node *t) /* ptr to node in a list  */
    ;
void *listAlloc(		   /* returns calloc() value */
		unsigned int n,	   /* number of items */
		unsigned int size) /* size of items */
    ;
void listNodeCheck(		      /* no return */
		   struct Node *node, /* Node to be checked */
		   const char *text)	      /* descriptive text about Node */
    ;
void listNodeDelete(		       /* no return */
		    struct Node *node, /* Node to be deleted */
		    const char *text)       /* descriptive text about node */
    ;
void listDeleteListList(struct Node *list) /* list of lists of items */
    ;
void listDeleteList(struct Node *list) /* list of items */
    ;
int listDeleteNext(struct Node *t) /* ptr to node in a list */
    ;
void listFree(		   /* no return */
	      void *p,	   /* pointer to segment to be freed */
	      const char *text) /* descriptive text about segment */
    ;

int mathZernike(	     /* returns #coeffs in u[] & ui[]    */
		double r,    /* radius in range 0-->1.0          */
		double cost, /* cosine of theta                  */
		double sint, /* sine   of theta                  */
		int nord[],  /* order limits                     */
		double u[],  /* Zernike coefficients return here */
		int ui[])    /* Indicies of coeffs return here   */
    ;

int mathSyminv2(	    /* returns TRUE if matrix is singular    */
		double a[], /* input a[n][n], overwritten by inverse */
		int n,	    /* order of 2-D matrix                   */
		int dim2)   /* size of fastest changing axis         */
    ;

#ifdef __cplusplus
};
#endif /* endif cplusplus */
#endif /* RAY_H */
