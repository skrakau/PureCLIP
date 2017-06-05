# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Generate Release builds by default
if ( NOT CMAKE_BUILD_TYPE )
    message ( STATUS "Setting build type to Release" )
    set ( CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type" FORCE )
else()
    message ( STATUS "User selected build type: ${CMAKE_BUILD_TYPE}" )
endif()

# Name of project and that it is C++ only.
project (PureCLIP CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

LIST ( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ )
find_package ( OpenMP REQUIRED )
include ( SeqAn )
include ( Boost )
find_package (GSL REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${GSL_INCLUDE_DIRS})
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


# PureCLIP static build for Linux
if ( CMAKE_PURECLIP_LINUX_STATIC )
    message ( STATUS "Configuring Linux static build" )
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()

# Add executable and link against SeqAn dependencies.
add_executable (pureclip
                    pureclip.cpp
                    util.h
                    call_sites.h
                    parse_alignments.h
                    prepro_mle.h
                    hmm_1.h
                    density_functions.h)

add_executable (winextract winextract.cpp)

target_link_libraries (pureclip ${Boost_LIBRARIES} ${SEQAN_LIBRARIES} ${GSL_LIBRARIES})
target_link_libraries (winextract ${Boost_LIBRARIES} ${SEQAN_LIBRARIES} ${GSL_LIBRARIES})
