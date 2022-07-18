cmake_minimum_required(VERSION 3.22)
project(optima C)

set(CMAKE_C_STANDARD 11)



include_directories(src)

add_executable(optima
        src/CPUCode/boundaries.c
        src/CPUCode/boundaries.h
        src/CPUCode/lb.c
        src/CPUCode/lb.h
        src/CPUCode/unsteady.c)

target_link_libraries(${PROJECT_NAME} m)
