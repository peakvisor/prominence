cmake_minimum_required(VERSION 3.10.0)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall -ffast-math -Wno-missing-braces")
set(CMAKE_CXX_STANDARD 17)

project(pcal)
project(fileprep)
project(pp)
project(geomorph)
project(pondem)
project(tests)

add_subdirectory(lib)

add_executable(pcal prominence_calculation/prominence_calculation.cpp prominence_calculation/isolation.cpp prominence_calculation/qtree.cpp prominence_calculation/qtree_node.cpp peaks/peak.cpp)
target_include_directories(pcal PUBLIC lib dem prominence_calculation peaks)
target_link_libraries(pcal index)

add_executable(fileprep dem/dem_processing.cpp)
target_link_libraries(fileprep index)
target_include_directories(fileprep PUBLIC lib dem)

add_executable(pp peaks/peaks_processing.cpp peaks/peak.cpp)
target_include_directories(pp PUBLIC lib peaks dem)
target_link_libraries(pp index)

add_executable(geomorph geomorphology/geomorphology.cpp peaks/peak.cpp)
target_include_directories(geomorph PUBLIC geomorphology lib peaks dem)
target_link_libraries(geomorph index)

add_executable(pondem peaks_on_dem_processing/peaks_on_dem_processing.cpp)
target_include_directories(pondem PUBLIC peaks_on_dem_processing dem)

add_executable(tests tests/tests.cpp prominence_calculation/isolation.cpp prominence_calculation/qtree.cpp prominence_calculation/qtree_node.cpp peaks/peak.cpp)
target_include_directories(tests PUBLIC lib lib dem peaks prominence_calculation geomorphology tests)
target_link_libraries(tests index)