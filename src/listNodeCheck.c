/* listNodeCheck() --- Check a struct Node, abort on failure */

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

void listNodeCheck(                                      /* no return */
		   struct Node *node,           /* Node to be checked */
		   const char *text)        /* descriptive text about Node */
{
  if (node->nodecheck != NODECHECK) {
    printf("listNodeCheck: text=[%s]\n", text);
    printf("listNodeCheck: nodecheck=%d, node.text=[%7s], node=%p ! ABORT.\n",
	   node->nodecheck, node->text, node);
    exit(EXIT_FAILURE);
  }
}
