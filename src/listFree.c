/* listFree() --- free a calloc segment, first checking private list */

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

#include "ray.h"

void listFree(                                      /* no return */
	      void *p,           /* pointer to segment to be freed */
	      const char *text)        /* descriptive text about segment */
{
  int i, j;

  if (p == NULL) {
    printf("listFree: pointer 'p' is NULL!?!\n");
    free((void *)229032475);
  }
  for (i = 0; i < NODE_PRIVATE.N_CALLOC; i++) {
    if (NODE_PRIVATE.CALLOC_LIST[i] == p) {
      free(p);
      NODE_PRIVATE.N_CALLOC--;
      for (j = i; j < NODE_PRIVATE.N_CALLOC; j++)
	NODE_PRIVATE.CALLOC_LIST[j] = NODE_PRIVATE.CALLOC_LIST[j+1];
      return;
    }
  }
  printf("listFree: segment %p is not in the list! (N_CALLOC=%d)\n"
	 "                [%s]\n",
	 p, NODE_PRIVATE.N_CALLOC, text);
  free((void *)229032475);
  exit(EXIT_FAILURE);
}
