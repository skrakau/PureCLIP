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
set (LINK_TYPE "-dynamic")
if ( STATIC_BUILD )
    set (LINK_TYPE "-static")
    if(UNIX AND NOT APPLE)
        message ( STATUS "Configuring Linux static build" )
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    elseif (UNIX AND APPLE)
        set (LINK_TYPE "")
        message ( STATUS "Configuring Mac build with static libstdc++" )
        string (REPLACE "dylib" "a" GSL_LIBRARIES "${GSL_LIBRARIES}")
        SET(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
        SET(BUILD_SHARED_LIBRARIES OFF)
    endif()
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

# Installation
if ( PKG_BUILD )
    SET ( BINDIR "." )
    SET ( EXTRADIR "." )
else()
    SET ( BINDIR "bin" )
    SET ( EXTRADIR "share/pureclip" )
endif()
install (TARGETS pureclip winextract RUNTIME DESTINATION "${BINDIR}")
install (FILES "../LICENSE.md" DESTINATION "${EXTRADIR}")

# Packaging
if(UNIX AND NOT APPLE)
    SET(CPACK_SYSTEM_NAME "linux64${LINK_TYPE}")
elseif (UNIX AND APPLE)
    SET(CPACK_SYSTEM_NAME "mac64${LINK_TYPE}")
endif()
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "3")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_INSTALLED_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../util;util" "${CMAKE_CURRENT_SOURCE_DIR}/../pkg;.")
message(STATUS "Adding file ${CMAKE_CURRENT_SOURCE_DIR}/../pkg/README")
include ( CPack )
