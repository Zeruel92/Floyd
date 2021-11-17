cmake_minimum_required(VERSION 3.20)
project(Floyd C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -DMAX_ITERATIONS=100")
set(CMAKE_C_FLAGS_DEBUG "-D_DEBUG -O3 -Wall -Wextra -DMAX_ITERATIONS=1")
set(CMAKE_C_COMPILER mpiicc)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(floyd floyd.c)
target_link_libraries(floyd ${MPI_LIBRARIES})
