# Makefile.del -- Makefile for exported 'ray' (ray-tracing) library.
#
# This raytracing library has been constructed as a part of 
# NRAO's Green Bank Telescope project.
# This makefile also compiles the rayMain() program.
# The 'list' package of six C functions is included in this library.
#
# D.Wells, NRAO-CV
# 1995-06-dd: initial version?
# 1997-05-01: many name changes
# 1998-10-05: revised this Makefile for clean exportable version
# 1998-10-25: added listNodeCreate.o to OBJECTS
#
# -=-=-=-=-=-=-=- Begin Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Copyright (C) 1998 Associated Universities, Inc. Washington DC, USA.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Correspondence concerning GBT software should be addressed as follows:
# GBT Operations, National Radio Astronomy Observatory, P. O. Box 2,
# Green Bank, WV 24944-0002 USA
#
# -=-=-=-=-=-=-=- End Copyright Notice -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# Substitute the path and options for your ANSI-C compiler here:
#    (XO5 selects highest level of optimization)
# CFLAGS	= -xO5
# ANSICC	= /opt/SUNWspro/bin/cc ${CFLAGS}
#.c.o:
#	-${ANSICC} -c $< 
#
LIBRAY	= libray.a
OBJECTS	= rayTrace.o rayAddSurface.o rayPrtSystem.o rayGenerator.o \
		rayPrtBundles.o rayPrtSegments.o rayPltPS.o rayPltSystem.o \
		rayGetFoci.o rayPrtFoci.o rayGetPlanes.o rayPrtPlanes.o \
		listAppend.o listDeleteList.o listDeleteListList.o \
		listDeleteNext.o listInitialize.o listInsertAfter.o \
		listNodeCreate.o listNodeDelete.o listNodeCheck.o \
		listAlloc.o listFree.o mathZernike.o mathSyminv2.o
TESTS	= rayTestCase1.in rayTestCase2.in rayTestCase3.in rayTestCase4.in 
OUTS	= ${TESTS:.in=.out}
DIFS	= ${TESTS:.in=.dif}
#
# Default make target produces ray tracing library and rayMain executable,
#     and then executes four test case:
rayPackage:	${LIBRAY} rayMain ${DIFS}
#
# Dependencies:
${OBJECTS} rayMain.o:	ray.h mathVectorMatrix.h gbtOpticalConstants.h
#
# Compile ray tracing library (commands may be operating system dependent):
${LIBRAY}:	${OBJECTS} 
	rm -f $@
	ar cqv $@ ${OBJECTS}
#
# Link main program (options may be operating system dependent):
rayMain:	${LIBRAY} rayMain.o 
	cc        -o $@ rayMain.o -L. -lray -lm
#	${ANSICC} -o $@ rayMain.o -L. -lray -lm
#
# Execute main program on a *.in test case:
%.out:	%.in rayMain
	rm -f ${*}*.ps
	-rayMain $< >${@} 2>&1
	cat ${@}
#
# Difference test case output against reference file %.ref:
%.dif:	%.out
	diff ${*}.out ${*}.ref | tee ${@}
#
clean_rayPackage:	
	rm -f core *~ *.o a.out
	rm -f ${LIBRAY} rayMain ${OUTS} rayTestCase*.ps ${DIFS}
