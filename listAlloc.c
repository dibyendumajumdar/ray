/* listAlloc() --- calloc() with private list of pointers */

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

void *listAlloc(                          /* returns calloc() value */
		unsigned int n,           /* number of items */
		unsigned int size)        /* size of items */
{
  void *p;
  
  if ((p = calloc(n, size)) == NULL) {
    printf("listAlloc: failed to calloc() n*size!?!\n");
    exit(EXIT_FAILURE);
  }
  NODE_PRIVATE.CALLOC_LIST[NODE_PRIVATE.N_CALLOC] = p;
  NODE_PRIVATE.N_CALLOC++;
  if ((NODE_PRIVATE.N_CALLOC < 0) ||
      (NODE_PRIVATE.N_CALLOC >= NODE_MAX)) {
    printf("listAlloc: N_CALLOC = %d. ABORT. (NODE_MAX=%d)\n",
	   NODE_PRIVATE.N_CALLOC, NODE_MAX);
    exit(EXIT_FAILURE);
  }
  return(p);
}
