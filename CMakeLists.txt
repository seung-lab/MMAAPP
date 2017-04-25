cmake_minimum_required (VERSION 2.8.11 FATAL_ERROR)

project (MMAAPP)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(Qt5Core)
find_package(Qt5Concurrent)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -fPIC")

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#add_subdirectory(examples)
add_definitions(-std=c++11 -DDEBUG -g -ggdb)
#add_definitions(-std=c++11 -O2 -g -ggdb -DDEBUG -D_USE_MATH_DEFINES -march=native)
add_definitions(-D_USE_MATH_DEFINES)

include_directories(
    ${Qt5Core_INCLUDE_DIRS}
    ${Qt5Concurrent_INCLUDE_DIRS}
    )

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(MMAAPP
    src/Info.cpp
    src/Segments.cpp
    src/main.cpp)

target_link_libraries(MMAAPP
    ${Qt5Core_LIBRARIES}
    ${Qt5Concurrent_LIBRARIES}
    )
