
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
the package, and input files for two example problems.  

## Copyright/LICENSE

The source files of library 'ray' and program 'rayMain' are
being made available under a **GNU "copyleft"** license. See file
"GNU_GPL_2.txt" in the ray distribution kit.  Each of the source files
begins with a copyright notice: 

> "Copyright (C) 1995 Associated
Universities, Inc. Washington DC, USA. This program is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation.."

## Notes

The software was obtained from NRAO ftp site. 
The original author is [Donald C. Wells](https://www.cv.nrao.edu/~dwells/).

## Changes 

* Moved to CMake build
* Changed eccetricity parameter to `k`, equaivalent to `1-eps^2` in original code
* Removed unused vignette parameters on surfaces
