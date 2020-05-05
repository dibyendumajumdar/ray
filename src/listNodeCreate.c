/* listNodeCreate() --- create a struct Node, return pointer to it */

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

struct Node *listNodeCreate(                /* return ptr to listAlloc(Node) */
			    char text[])       /* descriptive text re Node   */
{
  struct Node *node;
  
  if ((node = (struct Node *) listAlloc(1, sizeof(struct Node))) == NULL) {
    printf("listNodeCreate: could not listAlloc(1,Node) !?!\n");
    exit(EXIT_FAILURE);
  }
  if (NODE_DEBUG) printf("listNodeCreate: %7s  node=%8x\n", 
			 text, (unsigned int)node);
  node->nodecheck = NODECHECK;
  strncpy(node->text, text, NODETEXT-1);
  node->text[NODETEXT] = '\0'; 
  node->item = NULL;
  return (node);
}
