cmake_minimum_required(VERSION 3.22)
project(optima C)

set(CMAKE_C_STANDARD 11)



include_directories(C-bgk)

add_executable(optima
        C-bgk/boundaries.c
        C-bgk/boundaries.h
        C-bgk/lb.c
        C-bgk/lb.h
        C-bgk/unsteady.c)

target_link_libraries(${PROJECT_NAME} m)
