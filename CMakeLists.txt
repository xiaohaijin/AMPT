cmake_minimum_required(VERSION 3.0)

project(ampt)

enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMake_Fortran_FLAGS} -ffree-line-length-0")

add_subdirectory(src)
