cmake_minimum_required(VERSION 3.12)
project(ray VERSION 0.0.1 LANGUAGES C)

set(HEADERS
        ray.h
        mathVectorMatrix.h
        gbtOpticalConstants.h)

set(SRCS
        rayTrace.c
        rayAddSurface.c
        rayPrtSystem.c
        rayGenerator.c
        rayPrtBundles.c
        rayPrtSegments.c
        rayPltPS.c
        rayPltSystem.c
        rayGetFoci.c
        rayPrtFoci.c
        rayGetPlanes.c
        rayPrtPlanes.c
        listAppend.c
        listDeleteList.c
        listDeleteListList.c
        listDeleteNext.c
        listInitialize.c
        listInsertAfter.c
        listNodeCreate.c
        listNodeDelete.c
        listNodeCheck.c
        listAlloc.c
        listFree.c
        mathZernike.c
        mathSyminv2.c
        )

message("SOURCE dir is ${ray_SOURCE_DIR}")

if (WIN32)
    # disable warnings about C string functions
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif ()

if (NOT MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra -Wshadow -Wsign-compare -Wundef -Wwrite-strings \
        -Wredundant-decls -Wdisabled-optimization -Wdouble-promotion \
    	-Wlogical-op -Wno-aggressive-loop-optimizations \
        -Wmissing-prototypes -Wnested-externs \
        -Wstrict-prototypes -Werror=return-type -Wall -Wno-unused-parameter") # -Wc++-compat
endif ()

include(GNUInstallDirs)

set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

add_library(ray STATIC
        ${HEADERS}
        ${SRCS})

add_executable(raymain rayMain.c)
target_link_libraries(raymain ray)
