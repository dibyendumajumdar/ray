
		    The 'ray' Ray Tracing Package

			      Don Wells
			  <dwells@nrao.edu>

		      August 1995, October 1998


The "ray" package and program `rayMain` are capable of tracing sets of
rays through systems of rotationally-symmetric aspheric optical
elements.  The starting sets of rays can represent either plane or
spherical wavefronts, with feedhorn tapering. The optical elements can
be de-centered and/or tilted conic sections (spheres, ellipsoids,
paraboloids, hyperboloids) with additional superimposed
radially-symmetric aspheric terms, and they can be mirrors as well as
refracting surfaces.

This package includes full source for 'ray', plus a memo describing
the package, and input files for two example problems.  The kit
includes the `Makefile` which will build 'ray' and 'rayMain' and will
execute 'rayMain' on the test problems.

On a Sun, it should be possible to just do a 'make' in the directory.
On other Unix machines, it should only be necessary to change the CC
definition. On PCs, it may be necessary to transliterate the Makefile,
or to hand-execute equivalent operations. [NOT YET TESTED]

## Copyright/LICENSE

The source files of library 'ray' and program 'rayMain' are
being made available under a **GNU "copyleft"** license. See file
"GNU_GPL_2.txt" in the ray distribution kit.  Each of the source files
begins with a copyright notice: 

> "Copyright (C) 1995 Associated
Universities, Inc. Washington DC, USA. This program is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation.."

## NOTE ADDED 1998-10-25

C function `rayGetPlanes()` could cause
segmentation faults because it was erroneously calling function `free()`
with pointers to addresses which had not been allocated by `malloc()`.
The calls to `free()` have been deleted. In the course of searching for
this bug, several new functions were added to the list processing
package to create, delete and check the validity of list nodes. Also
function `rayPltPS()` was changed to permit overwriting previously
existing output files. Several other erroneous `free()` operations were
found, and then a systematic search found all memory leaks. -Don Wells

## Additional Notes

The software was obtained from NRAO ftp site. 
The original author is [Donald C. Wells](https://www.cv.nrao.edu/~dwells/).
