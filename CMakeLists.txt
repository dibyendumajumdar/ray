cmake_minimum_required(VERSION 3.12)
project(ray VERSION 0.0.1 LANGUAGES C)

set(HEADERS
        src/ray.h
        src/mathVectorMatrix.h
        src/gbtOpticalConstants.h)

set(SRCS
        src/rayTrace.c
        src/rayAddSurface.c
        src/rayPrtSystem.c
        src/rayGenerator.c
        src/rayPrtBundles.c
        src/rayPrtSegments.c
        src/rayPltPS.c
        src/rayPltSystem.c
        src/rayGetFoci.c
        src/rayPrtFoci.c
        src/rayGetPlanes.c
        src/rayPrtPlanes.c
        src/listAppend.c
        src/listDeleteList.c
        src/listDeleteListList.c
        src/listDeleteNext.c
        src/listInitialize.c
        src/listInsertAfter.c
        src/listNodeCreate.c
        src/listNodeDelete.c
        src/listNodeCheck.c
        src/listAlloc.c
        src/listFree.c
        src/mathZernike.c
        src/mathSyminv2.c
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

# Misc setup
if (NOT WIN32)
    set(EXTRA_LIBRARIES m)
endif ()

include_directories("${RaviCompiler_SOURCE_DIR}/src")
add_library(ray STATIC
        ${HEADERS}
        ${SRCS})
target_link_libraries(ray ${EXTRA_LIBRARIES})

add_executable(raymain src/rayMain.c)
target_link_libraries(raymain ray)
