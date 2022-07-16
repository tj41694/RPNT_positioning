# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-std=c++11)

project ("RPNT_cross_plat")

# 包含子项目。
add_subdirectory ("GNSSNav")
