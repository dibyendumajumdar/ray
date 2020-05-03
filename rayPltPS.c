/* rayPltPostscript() -- produce encapsulated Postscript from Segment list
   1998-06-11: decided to return non-zero on file error, rather than exit()
   1998-10-25: recover Hue_Table segments
*/

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1998 Associated Universities, Inc. Washington DC, USA.

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

#include <time.h>
#include "ray.h"

#define C2P(p) ((72.0/2.54)*(p))
#define S2P(p,w)  ((width_cm/2.54*72.0)*((p)/width)+(w/2.0))
#define XY2P(V,i,j) \
                (S2P((V[i])-to_point[i],C2P(width_cm))), \
                (S2P((V[j])-to_point[j],C2P(height_cm)))

struct HueTable {
    int ColorCode;  /* integer code number from Segment file */
    double hue;     /* random number assigned for Hue code   */
};

static struct Node *Hue_Table = NULL;

int rayPltPS(                       /* returns non-zero on error     */
	     struct Node *segments, /* list of line segments         */
	     double height_cm,      /* height of PS in centimeters   */
	     double width_cm,       /* width  of PS in centimeters   */
	     double width,          /* width in units of System      */
	     double to_point[],     /* XYZ of point in center of PS  */
	     char plt_mode[],       /* "Orthographic"|"Perspective"  */
	     char psname[])         /* Postscript output file        */
{
    int i, prev_color, active_path, jx, jy;
    const int eps=1, base_error_code=1100;
    double prev_T[3], diff;
    struct Node *s, *h;
    struct Segment *Seg;
    struct HueTable *hue;
    FILE *fp;
    time_t now;
    struct tm *now_tm;
    char header[FORMMAX], temp[FORMMAX];
    
    if ((fp = fopen(psname, "r")) != NULL) {   /* test if psname[] exists */
	printf("# rayPltPS WARNING: "
	       "Output file <%s> already exists, will overwrite.\n", 
	       psname);
	/* return(base_error_code+1); */
    }
    if ((fp = fopen(psname, "w")) == NULL) {   
	printf("rayPltPostscript: output file <%s> could not be created!?!\n",
	       psname);
	return(base_error_code+2);
    }
    if (segments == NULL) {
	printf("rayPltPS: pointer to <segments> is NULL!?!\n");
	return(base_error_code+3);
    }
    
    fprintf(fp, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(fp, "%%%%BoundingBox: %d %d %d %d\n",
	    0, 0, (int)ceil(C2P(width_cm)), (int)ceil(C2P(height_cm)));
    fprintf(fp, "%%%%BeginDocument: %s\n",
	    psname);
    fprintf(fp, "%%%%Creator: %s@%s\n",
	    (char *)getenv("LOGNAME"), (char *)getenv("DOMAIN"));
    time(&now);
    now_tm = gmtime(&now);
    fprintf(fp, "%%%%CreationDate: UTC %s",
	    asctime(now_tm));
    fprintf(fp, "%%%%Title: %s\n",
	    psname);
    fprintf(fp, "%%%%Pages: 1 1\n");
    fprintf(fp, "%%%%DocumentFonts: Courier\n");
    fprintf(fp, "%%%%EndComments\n");
    fprintf(fp, "%%%%BeginProcSet: graphics.pro\n\n");
    fprintf(fp, "%%%%EndProcSet\n");
    fprintf(fp, "%%%%EndProlog\n\n");
    
    fprintf(fp, "newpath 0 0 moveto %6.1f 0 rlineto 0 %6.1f rlineto\n",
	    C2P(width_cm), C2P(height_cm));
    fprintf(fp, "    %6.1f 0 rlineto 0 %6.1f rlineto closepath clip\n\n",
	    C2P(-width_cm), C2P(-height_cm));
    
    if (Hue_Table == NULL)
	Hue_Table = listInitialize("ColorCode-to-HSB"); /* temporary list */
    
    for (i = 0; i < 3; i++) prev_T[i] = -9999999.9;
    prev_color = -999999;
    active_path = FALSE;
    s = segments;
    while (s->next->next != s->next) {         /* loop over list of segments */
	s = s->next;
	Seg = (struct Segment *) s->item;
	
	if (Seg->ColorCode%1000 != prev_color) {
	    if (active_path) {
		fprintf(fp, "stroke\n%d setlinewidth\n",
			(Seg->ColorCode)/1000);
		active_path = FALSE;
	    }
	    h = Hue_Table;
	    hue = NULL;
	    while (h->next->next != h->next) { /* search for ColorCode */
		h = h->next;
		if (((struct HueTable *)h->item)->ColorCode
		    == Seg->ColorCode%1000) {
		    hue = (struct HueTable *)h->item;
		    break;
		}
	    }
	    if (hue == NULL) { /* we didn't find ColorCode, so create it */
		hue = (struct HueTable *) listAlloc(1, sizeof(struct HueTable));
		hue->ColorCode = Seg->ColorCode%1000;
		for (i = 1; i < 17; i++) {
		    hue->hue = (double)rand()/(double)RAND_MAX;
		    if ((((hue->hue > 0.30) && (hue->hue < 0.40)) && 1) ||
			(((hue->hue > 0.60) && (hue->hue < 0.99)) && 1)) break;
		} 
		/* fprintf(stdout, "Code=%3d  Hue=%4.2f  (i=%2d)\n",
			hue->ColorCode, hue->hue, i); */
		listAppend(hue, Hue_Table);
	    }
	    fprintf(fp, "%5.3f 1 1 sethsbcolor\n", hue->hue);
	    prev_color = Seg->ColorCode;
	    fprintf(fp, "newpath\n");
	    active_path = TRUE;
	    for (i = 0; i < 3; i++) prev_T[i] = -9999999.9;
	}
	for (i = 0, diff = 0.0; i < 3; i++) diff += fabs(Seg->T1[i] - prev_T[i]);
	jx = 2; jy = 0;
	jx = 0; jy = 1;
	if (diff != 0.0) fprintf(fp, "%8.1f %8.1f moveto\n",
				 XY2P(Seg->T1, jx, jy));
	fprintf(fp, "%8.1f %8.1f lineto\n", XY2P(Seg->T2, jx, jy));
	for (i = 0; i < 3; i++) prev_T[i] = Seg->T2[i];
    }
    if (active_path) fprintf(fp, "stroke\n");
    
    fprintf(fp, "0 0 0 sethsbcolor\n");
    
    fprintf(fp, "newpath 0 0 moveto %6.1f 0 rlineto 0 %6.1f rlineto\n",
	    C2P(width_cm)-eps, C2P(height_cm)-eps);
    fprintf(fp, "    %6.1f 0 rlineto 0 %6.1f rlineto closepath\n",
	    C2P(-width_cm)+eps, C2P(-height_cm)+eps);
    fprintf(fp, "    0.5 setlinewidth stroke\n\n");
    
    strcpy(header, (char *)segments->item);
    sprintf(temp, "; width=%.2f at [%.1f,%.1f]",
	    width, to_point[0], to_point[1]);
    strcat(header, temp);
    fprintf(fp, "/Helvetica-Narrow findfont 9 scalefont setfont\n");
    fprintf(fp, "%8.1f %8.1f moveto (%s) show\n",
	    0.0+C2P(0.1), C2P(height_cm)-C2P(0.35), header);
    
    fprintf(fp, "showpage\n");
    fprintf(fp, "%%%%EOF\n");
    fprintf(fp, "%%%%EndDocument\n");
    fclose(fp);

    return(0);
}

void rayPltPSDeleteHues()       /* no args, no return value */
{
  int i;
  struct Node *h;

  i = 0;
  h = Hue_Table;
  while (h->next->next != h->next) { /* count ColorCodes */
    h = h->next;
    i++;
  }
  /* printf("# rayPltPSDeleteHues: %d ColorCode nodes in list\n", i); */
  listDeleteList(Hue_Table); /* 1998-10-25: recover this memory! */
}
