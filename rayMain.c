/* rayMain() --- main program which calls functions of |ray| library
   D.Wells, NRAO-Cv, July95. */

/* -=-=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-

   Copyright (C) 1995 Associated Universities, Inc. Washington DC, USA.

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

#define LONG 300
#define SHORT 30

int main(int argc, char *argv[])
{
	char str[LONG], line[LONG], cmd[SHORT], c2[SHORT], wave_type[SHORT], assign_by[SHORT], name[LONG],
	    vignette[SHORT], temp[LONG], filename[SHORT], form[LONG], print_mode[LONG], psname[SHORT], pltmode[LONG],
	    set_name[LONG];
	char *p, *p_set_name = NULL;
	double xyz[3], XYZ[3], radius, step, angle, db, c, e, a2, a4, mu1, v[3], V[3], tol = 0.00001, hcm, wcm, wps,
										       to[3];
	int nv, nve, steps, num, code1, code2, d = 4, i, j, k = 65, axis_mask;
	enum VignetteType vign_type;
	enum ColorType a_b;
	struct Node *BundleSet, *System, *Segments, *Temp_List, *Foci, *Planes;
	FILE *fp;

	if (argc == 1)
		strcpy(filename, "rayTestCase4.in"); /* default to stdin eventually */
	else
		strcpy(filename, *++argv);
	if ((fp = fopen(filename, "r")) == NULL) {
		printf("rayMain: can't open <%s>\n", filename);
		exit(EXIT_FAILURE);
	}
	strcpy(line, "");
	NODE_PRIVATE.N_CALLOC = 0;
	BundleSet = listInitialize("");
	System = listInitialize("[empty list of surfaces]");
	Segments = listInitialize("");
	Foci = listInitialize("[empty list of foci]");
	Planes = listInitialize("[empty list of planes]");

	while (fgets(str, LONG, fp) != NULL) {
		i = strlen(line);
		strncat(line, str, (j = strlen(str) - 1));
		line[i + j] = '\0';
		if (line[i + j - 1] == '\\') { /* is it a continuation line? */
			line[i + j - 1] = '\0';
			strcat(line, " ");
			continue;
		}
		if (line[0] == '#') { /* is it a comment line? */
			printf("%s\n", line);

		} else { /* edit the line before scanning it: */
			while ((p = strpbrk(line, "\t,(){}[]\"&@")) != NULL) {
				i = (p - &line[0]); /* punctuation-->blanks */
				strncpy(temp, line, i);
				temp[i] = '\0';
				strcat(temp, " ");
				strcat(temp, ++p);
				strcpy(line, temp);
			}
			while ((p = strstr(line, "  ")) != NULL) {
				i = (p - &line[0]); /* twoblanks-->oneblank */
				strncpy(temp, line, i);
				temp[i] = '\0';
				strcat(temp, ++p);
				strcpy(line, temp);
			}
			while (line[0] == ' ') { /* delete leading blanks */
				strcpy(temp, &line[1]);
				strcpy(line, temp);
			}
			while (line[(i = strlen(line)) - 1] == ' ') {
				strncpy(temp, line, --i); /* delete trailing blanks */
				temp[i] = '\0';
				strcpy(line, temp);
			}
			if (strcmp(line, "") == 0) { /* is it a blank line? */
						     /* skip blank lines */

			} else { /* it is a command line: */
				if (sscanf(line, "%s ", cmd) != 1) {
					printf("line=<%s>\nLine is not empty, but no command code!?!\n", line);
					exit(EXIT_FAILURE);
				}
				printf("# Command=<"); /* print the command line */
				i = strlen(line);
				for (j = 0; j < (i - k); j += k) {
					strncpy(temp, &line[j], k);
					temp[k] = '\0';
					printf("%s|\n#         |", temp);
				}
				printf("%s>\n", &line[j]);

				if (strcmp(cmd, "Quit") == 0) {
					break;

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "System") == 0) {
					if ((nv = sscanf(line, "%s %s", c2, name)) != 2) {
						printf("Expected 2 values, got nv=%d -->ABORT\n", nv);
						exit(EXIT_FAILURE);
					}
					listDeleteList(System);
					System = listInitialize(name);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "Digits") == 0) {
					if ((nv = sscanf(line, "%s %d %lf", c2, &d, &tol)) != 3) {
						printf("Expected 3 values, got nv=%d -->ABORT\n", nv);
						exit(EXIT_FAILURE);
					}

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayAddSurface") == 0) {
					strcpy(form, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf");
					strcat(form, " %s %lf %lf %lf %lf %lf %lf %lf");
					if ((nv = sscanf(line, form, c2, name, &c, &e, &a2, &a4, &mu1, &xyz[0], &xyz[1],
							 &xyz[2], &XYZ[0], &XYZ[1], &XYZ[2], vignette, &v[0], &v[1],
							 &v[2], &V[0], &V[1], &V[2], &radius)) != 21) {
						printf("Expected 21 values, got nv=%d -->ABORT\n", nv);
						exit(EXIT_FAILURE);
					}
					if (strcmp(vignette, "cylinder") == 0)
						vign_type = VIGN_CYLINDER;
					else if (strcmp(vignette, "cone") == 0)
						vign_type = VIGN_CONE;
					else {
						printf("vignette_type_code=<%s> not recognized -->ABORT!\n", vignette);
						exit(EXIT_FAILURE);
					};
					rayAddSurface(System, name, c, e, a2, a4, mu1, xyz, XYZ, vign_type, v, V,
						      radius);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPrtSystem") == 0) {
					rayPrtSystem(System, d);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "SetName") == 0) {
					if ((nv = sscanf(line, "%s %s", c2, set_name)) != (nve = 2)) {
						printf("Expected %d values, got nv=%d -->ABORT\n", nve, nv);
						exit(EXIT_FAILURE);
					}
					p_set_name = &set_name[0];

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayGenerator") == 0) {
					strcpy(form, "%s %s %lf %lf %lf %lf %lf %lf %lf");
					strcat(form, "%lf %d %d %d %lf %lf %d %d %s");
					if ((nv = sscanf(line, form, c2, wave_type, &xyz[0], &xyz[1], &xyz[2], &XYZ[0],
							 &XYZ[1], &XYZ[2], &radius, &step, &steps, &axis_mask, &num,
							 &angle, &db, &code1, &code2, assign_by)) != (nve = 18)) {
						printf("Expected %d values, got nv=%d -->ABORT\n", nve, nv);
						exit(EXIT_FAILURE);
					}
					if (strcmp(assign_by, "bundle") == 0)
						a_b = COLOR_BUNDLE;
					else if (strcmp(assign_by, "ray") == 0)
						a_b = COLOR_RAY;
					else {
						printf("assign_by_type_code=<%s> not recognized -->ABORT!\n",
						       assign_by);
						exit(EXIT_FAILURE);
					};
					rayGenerator(BundleSet, p_set_name, wave_type, xyz, XYZ, radius, step, steps,
						     axis_mask, num, angle, db, code1, code2, a_b);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPrtBundles") == 0) {
					rayPrtBundles(BundleSet, d);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayTrace") == 0) {
					if ((nv = sscanf(line, "%s", c2)) != 1) {
						printf("Expected 1 values, got nv=%d -->ABORT\n", nv);
						exit(EXIT_FAILURE);
					}
					Temp_List = rayTrace(BundleSet, System, tol, Segments);
					listDeleteListList(BundleSet);
					BundleSet = Temp_List;

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPrtSegments") == 0) {
					rayPrtSegments(Segments, d);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayGetFoci") == 0) {
					if ((nv = sscanf(line, "%s", c2)) != 1) {
						printf("Expected 1 value, got nv=%d -->ABORT\n", nv);
						exit(EXIT_FAILURE);
					}
					listDeleteList(Foci);
					Foci = rayGetFoci(BundleSet, tol, Segments);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPrtFoci") == 0) {
					rayPrtFoci(Foci, d);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayGetPlanes") == 0) {
					if ((nv = sscanf(line, "%s %lf %lf %lf %s", c2, &xyz[0], &xyz[1], &radius,
							 print_mode)) != 5) {
						printf("Expected 5 values, got nv=%d -->ABORT\n", nv);
						exit(EXIT_FAILURE);
					}
					listDeleteList(Planes);
					Planes = rayGetPlanes(BundleSet, xyz, radius, print_mode);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPrtPlanes") == 0) {
					rayPrtPlanes(Planes, d);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPltPS") == 0) {
					if ((nv = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %s %s", c2, &hcm, &wcm, &wps,
							 &to[0], &to[1], &to[2], pltmode, psname)) != (nve = 9)) {
						printf("Expected %d values, got nv=%d -->ABORT\n", nve, nv);
						exit(EXIT_FAILURE);
					}
					rayPltPS(Segments, hcm, wcm, wps, to, pltmode, psname);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else if (strcmp(cmd, "rayPltSystem") == 0) {
					if ((nv = sscanf(line, "%s", c2)) != (nve = 1)) {
						printf("Expected %d values, got nv=%d -->ABORT\n", nve, nv);
						exit(EXIT_FAILURE);
					}
					rayPltSystem(System, Segments);

					/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- */
				} else {
					printf("Command=<%s> not recognized -->ABORT!\n", line);
					exit(EXIT_FAILURE);
				}
			}
		}
		strcpy(line, "");
	};
	listDeleteList(Planes);
	listDeleteList(Foci);
	listDeleteListList(BundleSet);
	listDeleteList(System);
	listDeleteList(Segments);
	/* rayPltPSDeleteHues(); <---CAUSES SEGMENTATION FAULTS!?! */
	if (NODE_PRIVATE.N_CALLOC > 0) {
		printf("# rayMain: memory not fully recovered:\n"
		       "#             NODE_PRIVATE.N_CALLOC = %d\n",
		       NODE_PRIVATE.N_CALLOC);
		exit(EXIT_FAILURE);
	}
	exit(EXIT_SUCCESS);
}
