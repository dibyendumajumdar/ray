/* listDeleteListList() -- delete list of lists of nodes and malloc-ed items 
   D.Wells, NRAO-CV
   1995-xx-xx: original version
   1998-10-25: fixed nasty bug---listDeleteList() would delete a segment
                   and then listDeleteNext() would try to delete it again.
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

#include "ray.h"

void listDeleteListList (struct Node *list)    /* list of lists of items */
{
  if (list == NULL) {
      printf("listDeleteListList: pointer to list of lists is NULL!?!\n");
      exit(EXIT_FAILURE);
  }
  while (list->next->next != list->next) {  /* loop until no more lists */
    listDeleteList(list->next->item);       /* delete next list item   */
    list->next->item = listAlloc(1, sizeof(int)); /* allocate segment so
						     listDeleteNext()
						     can delete it! */
    if (! listDeleteNext(list)) {           /* delete next Node struct */
	printf("listDeleteListList: failed to delete next node!?!\n");
	exit(EXIT_FAILURE);
    }
  }
  listNodeDelete(list->next, "DLL.tail");   /* delete the tail Node     */
  listFree(list->item, "l.D.L.L: head text"); /* delete the head string[] */
  listNodeDelete(list, "l.D.L.L: head");         /* delete the head Node     */
}
