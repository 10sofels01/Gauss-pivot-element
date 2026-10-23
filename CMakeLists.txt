cmake_minimum_required(VERSION 3.16)
project(gauss_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall)

add_executable(a.out
    main.cpp
    matrix_io.cpp
    gauss_solver.cpp
)