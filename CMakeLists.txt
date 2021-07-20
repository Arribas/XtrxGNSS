# Copyright (C) 2021 C. Fernandez-Prades, CTTC cfernandez@cttc.es
# Copyright (C) 2021 Javier Arribas, CTTC jarribas@cttc.es
#


################################################################################
# Project setup
################################################################################
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(WARNING "In-tree build is bad practice. Try 'cd build && cmake ..'")
endif()
cmake_minimum_required(VERSION 3.5..3.19)
project(xtrx_gnss CXX C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

################################################################################
# XTRX_GNSS version information
################################################################################
set(THIS_IS_A_RELEASE OFF)   # only related to version name, no further implications.
if(NOT ${THIS_IS_A_RELEASE})
    # Get the current working branch
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the latest abbreviated commit hash of the working branch
    execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

set(VERSION_INFO_MAJOR_VERSION 0)
set(VERSION_INFO_API_COMPAT    9)
if(${THIS_IS_A_RELEASE})
    set(VERSION_INFO_MINOR_VERSION 0)
else()
    set(VERSION_INFO_MINOR_VERSION git-${GIT_BRANCH}-${GIT_COMMIT_HASH})
endif()

set(VERSION ${VERSION_INFO_MAJOR_VERSION}.${VERSION_INFO_API_COMPAT}.${VERSION_INFO_MINOR_VERSION})



################################################################################
# Environment setup
################################################################################
include(ExternalProject)
include(FeatureSummary)

# Detect Linux Distribution
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux|kFreeBSD|GNU")
    include(DetectLinuxDistro)
    if(CMAKE_CROSSCOMPILING)
        message(STATUS "Configuring XTRX_GNSS v${VERSION} to be cross-compiled on ${LINUX_DISTRIBUTION} ${LINUX_VER} (${CMAKE_HOST_SYSTEM_PROCESSOR}) for ${CMAKE_SYSTEM_PROCESSOR} ${ARCHITECTURE_STRING}")
    else()
        message(STATUS "Configuring XTRX_GNSS v${VERSION} to be built on GNU/Linux ${LINUX_DISTRIBUTION} ${LINUX_VER} ${ARCHITECTURE_STRING}")
    endif()
endif()

# Define extra build types and select Release by default to get optimization flags
include(XtrxGNSSBuildTypes)
# Available options:
#  - None: nothing set
#  - Debug: -O2 -g
#  - Release: -O3
#  - RelWithDebInfo: -O3 -g
#  - MinSizeRel: -Os
#  - Coverage: -Wall -pedantic -pthread -g -O0 -fprofile-arcs -ftest-coverage
#  - NoOptWithASM: -O0 -g -save-temps
#  - O2WithASM: -O2 -g -save-temps
#  - O3WithASM: -O3 -g -save-temps
#  - ASAN: -Wall -Wextra -g -O2 -fsanitize=address -fno-omit-frame-pointer
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to Release.")
else()
   message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}.")
endif()
xtrxgnss_check_build_type(${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "")

# Append -O2 optimization flag for Debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")

# allow 'large' files in 32 bit builds
if(UNIX)
    add_definitions(-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES)
endif()

# If this is an out-of-tree build, do not pollute the original source directory
if(${CMAKE_BINARY_DIR} MATCHES ${CMAKE_SOURCE_DIR})
    set(LOCAL_INSTALL_BASE_DIR ${CMAKE_SOURCE_DIR})
else()
    set(LOCAL_INSTALL_BASE_DIR ${CMAKE_BINARY_DIR})
endif()

# Determine if we are using make or ninja
if(CMAKE_MAKE_PROGRAM MATCHES "make")
    set(CMAKE_MAKE_PROGRAM_PRETTY_NAME "make")
endif()
if(CMAKE_MAKE_PROGRAM MATCHES "ninja")
    set(CMAKE_MAKE_PROGRAM_PRETTY_NAME "ninja")
endif()
if(NOT CMAKE_MAKE_PROGRAM_PRETTY_NAME)
    set(CMAKE_MAKE_PROGRAM_PRETTY_NAME "${CMAKE_MAKE_PROGRAM}")
endif()



################################################################################
# Set minimal C and C++ standards
################################################################################
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

# Visibility
# See https://gcc.gnu.org/wiki/Visibility
if(POLICY CMP0063)
    cmake_policy(SET CMP0063 NEW)
    set(CMAKE_CXX_VISIBILITY_PRESET hidden)
    set(CMAKE_C_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
else()
    if((CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND NOT WIN32)
        add_definitions(-fvisibility=hidden)
    endif()
endif()



################################################################################
# Minimum required versions
################################################################################
set(XTRX_GNSS_BOOST_MIN_VERSION "1.53")



################################################################################
# Versions to download and build (but not installed) if not found
################################################################################
set(XTRX_GNSS_GFLAGS_LOCAL_VERSION "2.2.2")
set(XTRX_GNSS_GLOG_LOCAL_VERSION "0.4.0")


################################################################################
# pkg-config - Helper tool used when compiling applications and libraries.
################################################################################
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH TRUE)
find_package(PkgConfig)


################################################################################
# Find the POSIX thread (pthread) libraries
################################################################################
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
if(CMAKE_CROSSCOMPILING)
    set(THREADS_PREFER_PTHREAD_FLAG FALSE)
else()
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
endif()
find_package(Threads REQUIRED)
set_package_properties(Threads PROPERTIES
    URL "https://computing.llnl.gov/tutorials/pthreads/"
    DESCRIPTION "Implements the POSIX Threads execution model"
    PURPOSE "Used to implement parallelism."
    TYPE REQUIRED
)



################################################################################
# Boost - https://www.boost.org
################################################################################
if(UNIX AND EXISTS "/usr/lib64")
    list(APPEND BOOST_LIBRARYDIR "/usr/lib64") # Fedora 64-bit fix
endif()
# Boost_ADDITIONAL_VERSIONS is only used internally by cmake to know the
# formation of newer versions. No need to increase, not used anymore since newer
# Boost provides its own CMake configuration files.
set(Boost_ADDITIONAL_VERSIONS
    "1.53.0" "1.53" "1.54.0" "1.54"
    "1.55.0" "1.55" "1.56.0" "1.56" "1.57.0" "1.57" "1.58.0" "1.58" "1.59.0" "1.59"
    "1.60.0" "1.60" "1.61.0" "1.61" "1.62.0" "1.62" "1.63.0" "1.63" "1.64.0" "1.64"
    "1.65.0" "1.65" "1.66.0" "1.66" "1.67.0" "1.67" "1.68.0" "1.68" "1.69.0" "1.69"
    "1.70.0" "1.70" "1.71.0" "1.71"
)
set(Boost_USE_MULTITHREAD ON)
set(Boost_USE_STATIC_LIBS OFF)
set(BOOST_COMPONENTS date_time system filesystem thread serialization chrono program_options)
find_package(Boost ${LIME_GNSS_BOOST_MIN_VERSION} COMPONENTS ${BOOST_COMPONENTS} REQUIRED)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Fatal error: Boost (version >=${LIME_GNSS_BOOST_MIN_VERSION}) required.")
endif()

set_package_properties(Boost PROPERTIES
    URL "https://www.boost.org"
    PURPOSE "Used widely across the source code."
    TYPE REQUIRED
)

if(CMAKE_VERSION VERSION_LESS 3.14)
    set(Boost_VERSION_STRING "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
endif()
if(POLICY CMP0093)
    cmake_policy(SET CMP0093 NEW)  # FindBoost reports Boost_VERSION in x.y.z format.
endif()

# Define Boost::headers target if CMake < 3.15
if(CMAKE_VERSION VERSION_LESS 3.15)
    if(NOT TARGET Boost::headers)
        add_library(Boost::headers INTERFACE IMPORTED)
        set_target_properties(Boost::headers PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR}
        )
    endif()
endif()

set_package_properties(Boost PROPERTIES
    DESCRIPTION "Portable C++ source libraries (found: v${Boost_VERSION_STRING})"
)



################################################################################
# gflags - https://github.com/gflags/gflags
################################################################################
set(LOCAL_GFLAGS FALSE)
find_package(GFLAGS)
set_package_properties(GFLAGS PROPERTIES
    PURPOSE "Used for commandline flags management."
    TYPE REQUIRED
)
if(NOT GFLAGS_FOUND)
    message(STATUS " gflags library has not been found.")
    message(STATUS " gflags v${LIME_GNSS_GFLAGS_LOCAL_VERSION} will be downloaded and built automatically")
    message(STATUS " when doing '${CMAKE_MAKE_PROGRAM_PRETTY_NAME}'.")
    set(GFLAGS_BUILD_COMMAND ${CMAKE_COMMAND}
        "--build" "${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}"
        "--config" $<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>$<$<CONFIG:NoOptWithASM>:Debug>$<$<CONFIG:Coverage>:Debug>$<$<CONFIG:O2WithASM>:RelWithDebInfo>$<$<CONFIG:O3WithASM>:RelWithDebInfo>$<$<CONFIG:ASAN>:Debug>
    )
    if(CMAKE_GENERATOR STREQUAL Xcode)
        set(GFLAGS_BUILD_COMMAND "xcodebuild" "-configuration" $<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>)
    endif()
    if(CMAKE_TOOLCHAIN_FILE)
        set(GFLAGS_TOOLCHAIN_FILE -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE})
    endif()

    set(GFLAGS_BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${CMAKE_STATIC_LIBRARY_SUFFIX})
    if((CMAKE_BUILD_TYPE STREQUAL Debug) OR (CMAKE_BUILD_TYPE STREQUAL NoOptWithASM) OR
        (CMAKE_BUILD_TYPE STREQUAL Coverage) OR (CMAKE_BUILD_TYPE STREQUAL ASAN))  # Workaround for Ninja generator
        set(GFLAGS_BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags_debug${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
    if((CMAKE_VERSION VERSION_GREATER 3.12.0) AND NOT (CMAKE_GENERATOR STREQUAL Xcode))
        set(PARALLEL_BUILD "--parallel 2")
    endif()
    ExternalProject_Add(gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}
        PREFIX ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}
        GIT_REPOSITORY git://github.com/gflags/gflags.git
        GIT_TAG v${LIME_GNSS_GFLAGS_LOCAL_VERSION}
        SOURCE_DIR ${CMAKE_BINARY_DIR}/thirdparty/gflags/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}
        BINARY_DIR ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}
        CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
            -DBUILD_STATIC_LIBS=ON
            -DBUILD_gflags_LIB=ON
            -DBUILD_gflags_nothreads_LIB=ON
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            ${GFLAGS_TOOLCHAIN_FILE}
            -DCMAKE_BUILD_TYPE=$<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>$<$<CONFIG:NoOptWithASM>:Debug>$<$<CONFIG:Coverage>:Debug>$<$<CONFIG:O2WithASM>:RelWithDebInfo>$<$<CONFIG:O3WithASM>:RelWithDebInfo>$<$<CONFIG:ASAN>:Debug>
        BUILD_COMMAND ${GFLAGS_BUILD_COMMAND} ${PARALLEL_BUILD}
        BUILD_BYPRODUCTS ${GFLAGS_BUILD_BYPRODUCTS}
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        INSTALL_COMMAND ""
    )
    # Note: -DBUILD_gflags_nothreads_LIB=ON is required as a workaround to a bug in gflags 2.2.2. This is fixed in gflags master branch

    set(GFlags_INCLUDE_DIRS
        ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/include CACHE PATH "Local Gflags headers"
    )

    if(NOT TARGET Gflags::gflags)
        file(MAKE_DIRECTORY ${GFlags_INCLUDE_DIRS})
        add_library(Gflags::gflags STATIC IMPORTED)
        add_dependencies(Gflags::gflags gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION})
        set_target_properties(Gflags::gflags PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
            IMPORTED_CONFIGURATIONS "None;Debug;Release;RelWithDebInfo;MinSizeRel"
            MAP_IMPORTED_CONFIG_NOOPTWITHASM Debug
            MAP_IMPORTED_CONFIG_COVERAGE Debug
            MAP_IMPORTED_CONFIG_O2WITHASM RelWithDebInfo
            MAP_IMPORTED_CONFIG_O3WITHASM RelWithDebInfo
            MAP_IMPORTED_CONFIG_ASAN Debug
            IMPORTED_LOCATION_NONE ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${CMAKE_STATIC_LIBRARY_SUFFIX}
            IMPORTED_LOCATION_DEBUG ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags_debug${CMAKE_STATIC_LIBRARY_SUFFIX}
            IMPORTED_LOCATION_RELEASE ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${CMAKE_STATIC_LIBRARY_SUFFIX}
            IMPORTED_LOCATION_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${CMAKE_STATIC_LIBRARY_SUFFIX}
            IMPORTED_LOCATION_MINSIZEREL ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${CMAKE_STATIC_LIBRARY_SUFFIX}
            INTERFACE_INCLUDE_DIRECTORIES ${GFlags_INCLUDE_DIRS}
            INTERFACE_LINK_LIBRARIES ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gflags$<$<CONFIG:Debug>:_debug>${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
        if((CMAKE_GENERATOR STREQUAL Xcode) OR MSVC)
            if(MSVC)
                set(MSVC_POSTFIX _static)
            endif()
            set_target_properties(Gflags::gflags PROPERTIES
                IMPORTED_LOCATION_DEBUG ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/Debug/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${MSVC_POSTFIX}_debug${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_RELEASE ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/Release/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${MSVC_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/RelWithDebInfo/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${MSVC_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_MINSIZEREL ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/MinSizeRel/${CMAKE_FIND_LIBRARY_PREFIXES}gflags${MSVC_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}
                INTERFACE_LINK_LIBRARIES ${CMAKE_BINARY_DIR}/gflags-${LIME_GNSS_GFLAGS_LOCAL_VERSION}/lib/$<$<CONFIG:Debug>:Debug/>$<$<CONFIG:Release>:Release/>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo/>$<$<CONFIG:MinSizeRel>:MinSizeRel/>${CMAKE_FIND_LIBRARY_PREFIXES}gflags${MSVC_POSTFIX}$<$<CONFIG:Debug>:_debug>${CMAKE_STATIC_LIBRARY_SUFFIX}
            )
        endif()
    endif()

    if(MSVC)
        target_link_libraries(Gflags::gflags INTERFACE shlwapi.lib)
    endif()

    set(LOCAL_GFLAGS TRUE CACHE STRING "GFlags downloaded and built automatically" FORCE)
    set_package_properties(GFLAGS PROPERTIES
        PURPOSE "Gflags v${LIME_GNSS_GFLAGS_LOCAL_VERSION} will be downloaded and built when doing '${CMAKE_MAKE_PROGRAM_PRETTY_NAME}'."
    )
endif()



################################################################################
# glog - https://github.com/google/glog
################################################################################
if(NOT ENABLE_OWN_GLOG AND NOT ${LOCAL_GFLAGS})
    find_package(GLOG)
endif()
set_package_properties(GLOG PROPERTIES
    PURPOSE "Used for runtime internal logging."
    TYPE REQUIRED
)
if(NOT GLOG_FOUND OR ${LOCAL_GFLAGS})
    message(STATUS " glog library has not been found")
    if(NOT GFLAGS_FOUND)
        message(STATUS " or it is likely not linked to gflags.")
    endif()
    message(STATUS " glog v${LIME_GNSS_GLOG_LOCAL_VERSION} will be downloaded and built automatically")
    message(STATUS " when doing '${CMAKE_MAKE_PROGRAM_PRETTY_NAME}'.")
    if(NOT ${LOCAL_GFLAGS})
        add_library(gflags-${GNSSSDR_GFLAGS_LOCAL_VERSION} UNKNOWN IMPORTED)
        set_property(TARGET gflags-${GNSSSDR_GFLAGS_LOCAL_VERSION} PROPERTY IMPORTED_LOCATION "${GFlags_LIBS}")
        string(REPLACE /include "" GFLAGS_PREFIX_PATH ${GFlags_INCLUDE_DIRS})
    else()
        set(GFLAGS_PREFIX_PATH ${CMAKE_BINARY_DIR}/gflags-${GNSSSDR_GFLAGS_LOCAL_VERSION})
    endif()
    set(TARGET_GFLAGS gflags-${GNSSSDR_GFLAGS_LOCAL_VERSION})
    set(GLOG_MAKE_PROGRAM ${CMAKE_COMMAND}
        "--build" "${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}"
        "--config" $<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>$<$<CONFIG:NoOptWithASM>:Debug>$<$<CONFIG:Coverage>:Debug>$<$<CONFIG:O2WithASM>:RelWithDebInfo>$<$<CONFIG:O3WithASM>:RelWithDebInfo>$<$<CONFIG:ASAN>:Debug>
    )
    if(CMAKE_GENERATOR STREQUAL Xcode)
        set(GLOG_MAKE_PROGRAM "xcodebuild" "-configuration"
            $<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>
        )
    endif()
    if(CMAKE_TOOLCHAIN_FILE)
        set(GLOG_TOOLCHAIN_FILE -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE})
    endif()

    if(CMAKE_VERSION VERSION_LESS 3.2)
        if(CMAKE_VERSION VERSION_LESS 3.0)
            set(GLOG_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM})
            set(GFLAGS_LIBRARIES_TO_LINK ${GFlags_LIBS})
            if(${LOCAL_GFLAGS})
                set(GFLAGS_LIBRARY_DIR_TO_LINK ${CMAKE_BINARY_DIR}/gflags-${GNSSSDR_GFLAGS_LOCAL_VERSION}/lib)
            else()
                set(GFLAGS_LIBRARY_DIR_TO_LINK ${GFlags_LIBRARY_DIRS})
            endif()
            if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
                set(GFLAGS_LIBRARIES_TO_LINK "${GFLAGS_LIBRARIES_TO_LINK} -lc++")
                set(GLOG_EXPORT_CXX_LIBRARIES "export CXXFLAGS=\"-stdlib=libc++\"")
            endif()
            if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                set(GLOG_EXPORT_C_COMPILER "export CC=clang")
                set(GLOG_EXPORT_CXX_COMPILER "export CXX=clang++")
            endif()
            file(WRITE ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/tmp/configure_with_gflags
"#!/bin/sh
export CPPFLAGS=-I${GFlags_INCLUDE_DIRS}
export LDFLAGS=-L${GFLAGS_LIBRARY_DIR_TO_LINK}
export LIBS=\"${GFLAGS_LIBRARIES_TO_LINK}\"
${GLOG_EXPORT_CXX_LIBRARIES}
${GLOG_EXPORT_C_COMPILER}
${GLOG_EXPORT_CXX_COMPILER}
cd ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/
aclocal
automake --add-missing
autoreconf -vfi
cd ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/configure --enable-shared=no"
            )

            file(COPY ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/tmp/configure_with_gflags
                DESTINATION ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
                    GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
            )

            set(GLOG_CONFIGURE ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/configure_with_gflags)

            # Ensure that aclocal and libtool are present
            if(${CMAKE_SYSTEM_NAME} MATCHES "Linux|kFreeBSD|GNU")
                if(EXISTS "/usr/bin/libtoolize")
                    if(EXISTS "/usr/bin/aclocal" OR
                        EXISTS "/usr/bin/aclocal-1.16" OR
                        EXISTS "/usr/bin/aclocal-1.15" OR
                        EXISTS "/usr/bin/aclocal-1.14" OR
                        EXISTS "/usr/bin/aclocal-1.13" OR
                        EXISTS "/usr/bin/aclocal-1.11" OR
                        EXISTS "/usr/bin/aclocal-1.10")
                        # Everything ok, we can move on
                    else()
                        message(" aclocal has not been found.")
                        message(" You can try to install it by typing:")
                        if(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
                            message(" sudo yum groupinstall 'Development Tools'")
                        elseif(${LINUX_DISTRIBUTION} MATCHES "openSUSE")
                            message(" sudo zypper install automake")
                        else()
                            message(" sudo apt-get install automake")
                        endif()
                        message(FATAL_ERROR "aclocal is required to build glog from source")
                    endif()
                else()
                    message(" libtool has not been found.")
                    message(" You can try to install it by typing:")
                    if(${LINUX_DISTRIBUTION} MATCHES "Fedora" OR ${LINUX_DISTRIBUTION} MATCHES "Red Hat")
                        message(" sudo yum groupinstall 'Development Tools'")
                    elseif(${LINUX_DISTRIBUTION} MATCHES "openSUSE")
                        message(" sudo zypper install libtoool")
                    else()
                        message(" sudo apt-get install libtool")
                    endif()
                    message(FATAL_ERROR "libtool is required to build glog from source")
                endif()
            endif()

            if(GLOG_MAKE_PROGRAM MATCHES "ninja")
                find_program(GLOG_MAKE_EXECUTABLE make
                    PATHS
                        /usr/bin
                        /usr/local/bin
                )
                if(NOT GLOG_MAKE_EXECUTABLE)
                    message(FATAL_ERROR "make is required to build Glog from source.")
                endif()
                set(GLOG_MAKE_PROGRAM ${GLOG_MAKE_EXECUTABLE})
            endif()
            ExternalProject_Add(glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                DEPENDS ${TARGET_GFLAGS}
                PREFIX ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                GIT_REPOSITORY https://github.com/google/glog/
                GIT_TAG v${LIME_GNSS_GLOG_LOCAL_VERSION}
                SOURCE_DIR ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                BINARY_DIR ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                CONFIGURE_COMMAND ${GLOG_CONFIGURE} --prefix=<INSTALL_DIR>
                BUILD_COMMAND "${GLOG_MAKE_PROGRAM}"
                UPDATE_COMMAND ""
                PATCH_COMMAND ""
                INSTALL_COMMAND ""
            )
            set(GLOG_LIBRARIES
                ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/.libs/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
            )
            set(GLOG_INCLUDE_DIRS
                ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/src
                ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/src
            )
        else()  # CMake > 3.0 but < 3.2
            ExternalProject_Add(glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                DEPENDS ${TARGET_GFLAGS}
                PREFIX ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                GIT_REPOSITORY https://github.com/google/glog/
                GIT_TAG v${LIME_GNSS_GLOG_LOCAL_VERSION}
                SOURCE_DIR ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                BINARY_DIR ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                    -DCMAKE_PREFIX_PATH=${GFLAGS_PREFIX_PATH}
                    ${GLOG_TOOLCHAIN_FILE}
                    -DCMAKE_BUILD_TYPE=$<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>$<$<CONFIG:NoOptWithASM>:Debug>$<$<CONFIG:Coverage>:Debug>$<$<CONFIG:O2WithASM>:RelWithDebInfo>$<$<CONFIG:O3WithASM>:RelWithDebInfo>$<$<CONFIG:ASAN>:Debug>
                BUILD_COMMAND ${GLOG_MAKE_PROGRAM}
                UPDATE_COMMAND ""
                PATCH_COMMAND ""
                INSTALL_COMMAND ""
            )
            set(GLOG_INCLUDE_DIRS
                ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/src
                ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
                ${GFlags_INCLUDE_DIRS}
            )
        endif()
    else()  # CMake > 3.2
        set(GLOG_BUILD_BYPRODUCTS
            ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
        )
        if((CMAKE_BUILD_TYPE STREQUAL Debug) OR (CMAKE_BUILD_TYPE STREQUAL NoOptWithASM) OR
            (CMAKE_BUILD_TYPE STREQUAL Coverage) OR (CMAKE_BUILD_TYPE STREQUAL ASAN))  # Workaround for Ninja generator
            set(GLOG_BUILD_BYPRODUCTS
                ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glogd${CMAKE_STATIC_LIBRARY_SUFFIX}
            )
        endif()
        if((CMAKE_VERSION VERSION_GREATER 3.12.0) AND NOT (CMAKE_GENERATOR STREQUAL Xcode))
            set(PARALLEL_BUILD "--parallel 2")
        endif()
        ExternalProject_Add(glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
            DEPENDS ${TARGET_GFLAGS}
            PREFIX ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
            GIT_REPOSITORY https://github.com/google/glog/
            GIT_TAG v${LIME_GNSS_GLOG_LOCAL_VERSION}
            SOURCE_DIR ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
            BINARY_DIR ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
            CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_PREFIX_PATH=${GFLAGS_PREFIX_PATH}
                ${GLOG_TOOLCHAIN_FILE}
                -DCMAKE_BUILD_TYPE=$<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo>$<$<CONFIG:MinSizeRel>:MinSizeRel>$<$<CONFIG:NoOptWithASM>:Debug>$<$<CONFIG:Coverage>:Debug>$<$<CONFIG:O2WithASM>:RelWithDebInfo>$<$<CONFIG:O3WithASM>:RelWithDebInfo>$<$<CONFIG:ASAN>:Debug>
            BUILD_COMMAND "${GLOG_MAKE_PROGRAM} ${PARALLEL_BUILD}"
            BUILD_BYPRODUCTS ${GLOG_BUILD_BYPRODUCTS}
            UPDATE_COMMAND ""
            PATCH_COMMAND ""
            INSTALL_COMMAND ""
        )
        set(GLOG_INCLUDE_DIRS
            ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/src
            ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}
            ${GFlags_INCLUDE_DIRS}
        )
    endif()

    add_dependencies(glog-${LIME_GNSS_GLOG_LOCAL_VERSION} Gflags::gflags)

    # Create Glog::glog target
    if(NOT TARGET Glog::glog)
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/thirdparty/glog/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/src)
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION})
        add_library(Glog::glog STATIC IMPORTED)
        add_dependencies(Glog::glog glog-${LIME_GNSS_GLOG_LOCAL_VERSION})
        if(CMAKE_VERSION VERSION_LESS 3.0)
            set_target_properties(Glog::glog PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                IMPORTED_LOCATION "${GLOG_LIBRARIES}"
                INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIRS}"
                INTERFACE_INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${GLOG_LIBRARIES}"
            )
        else()
            set_target_properties(Glog::glog PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                IMPORTED_CONFIGURATIONS "None;Debug;Release;RelWithDebInfo;MinSizeRel"
                MAP_IMPORTED_CONFIG_NOOPTWITHASM Debug
                MAP_IMPORTED_CONFIG_COVERAGE Debug
                MAP_IMPORTED_CONFIG_O2WITHASM RelWithDebInfo
                MAP_IMPORTED_CONFIG_O3WITHASM RelWithDebInfo
                MAP_IMPORTED_CONFIG_ASAN Debug
                IMPORTED_LOCATION_NONE ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_DEBUG ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glogd${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_RELEASE ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                IMPORTED_LOCATION_MINSIZEREL ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                INTERFACE_INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/${CMAKE_FIND_LIBRARY_PREFIXES}glog$<$<CONFIG:Debug>:d>${CMAKE_STATIC_LIBRARY_SUFFIX}
            )
            if((CMAKE_GENERATOR STREQUAL Xcode) OR MSVC)
                set_target_properties(Glog::glog PROPERTIES
                    IMPORTED_LOCATION_DEBUG ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/Debug/${CMAKE_FIND_LIBRARY_PREFIXES}glogd${CMAKE_STATIC_LIBRARY_SUFFIX}
                    IMPORTED_LOCATION_RELEASE ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/Release/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                    IMPORTED_LOCATION_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/RelWithDebInfo/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                    IMPORTED_LOCATION_MINSIZEREL ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/MinSizeRel/${CMAKE_FIND_LIBRARY_PREFIXES}glog${CMAKE_STATIC_LIBRARY_SUFFIX}
                    INTERFACE_LINK_LIBRARIES ${CMAKE_BINARY_DIR}/glog-${LIME_GNSS_GLOG_LOCAL_VERSION}/$<$<CONFIG:Debug>:Debug/>$<$<CONFIG:Release>:Release/>$<$<CONFIG:RelWithDebInfo>:RelWithDebInfo/>$<$<CONFIG:MinSizeRel>:MinSizeRel/>${CMAKE_FIND_LIBRARY_PREFIXES}glog$<$<CONFIG:Debug>:d>${CMAKE_STATIC_LIBRARY_SUFFIX}
                )
            endif()
        endif()
    endif()

    set(LOCAL_GLOG TRUE CACHE STRING "Glog downloaded and built automatically" FORCE)

    set_package_properties(GLOG PROPERTIES
        PURPOSE "Glog v${LIME_GNSS_GLOG_LOCAL_VERSION} will be downloaded and built when doing 'make'."
    )
endif()

if(NOT ENABLE_LOG)
    message(STATUS "Internal logging is not enabled")
    if(CMAKE_VERSION VERSION_GREATER 3.11.0)
        target_compile_definitions(Glog::glog INTERFACE -DGOOGLE_STRIP_LOG=1)
    else()
        set_property(TARGET Glog::glog APPEND PROPERTY
            INTERFACE_COMPILE_DEFINITIONS GOOGLE_STRIP_LOG=1
        )
    endif()
endif()


################################################################################
# Set compiler flags
################################################################################
set(CXX_WARNING_FLAGS -Wall -Wextra)
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND NOT WIN32)
    # Add warning flags
    # For "-Wall" see https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
    add_compile_options("$<$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,CXX>:${CXX_WARNING_FLAGS}>")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0")
        add_compile_options(-Wno-missing-field-initializers)
    endif()
    if(CMAKE_CROSSCOMPILING)
        add_compile_options(-Wno-psabi)
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options("$<$<STREQUAL:$<TARGET_PROPERTY:LINKER_LANGUAGE>,CXX>:${CXX_WARNING_FLAGS}>")
endif()


################################################################################
# Create uninstall target
################################################################################
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    @ONLY
)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)


################################################################################
# Add subdirectories (in order of deps)
################################################################################
#add_subdirectory(src)
add_subdirectory(tools)
