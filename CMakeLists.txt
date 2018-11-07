cmake_minimum_required(VERSION 3.0)

project(ampt CXX Fortran C)

enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(ROOT)
if(ROOT_FOUND)
  message(STATUS "查找到了需要使用的ROOT库")
else(ROOT_FOUND)
  message(ERROR "未查找到ROOT库，请查看下是否未安装或未配置")
endif(ROOT_FOUND)

add_subdirectory(src)
