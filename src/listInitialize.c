/* listInitialize() --- create head & tail nodes with descriptive string */

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

struct Node *listInitialize(                  /* return ptr to new list      */
			    char *name)       /* descriptive string for list */
{
  struct Node *head, *tail;

  head = listNodeCreate("head");
  if ((head->item = (char *) listAlloc(strlen(name)+1, sizeof(char))) == NULL) {
      printf("listInitialize: could not malloc head string!?!\n");
      exit(EXIT_FAILURE);
  }
  strcpy(head->item, name);    /* head nodes contain descriptive strings */
  tail = listNodeCreate("tail");
  head->next = tail;           /* initialize head node pointing to tail node */
  tail->item = NULL;           /* tail nodes have no data */
  tail->next = tail;           /* and point to themselves */
  return (head);
}
