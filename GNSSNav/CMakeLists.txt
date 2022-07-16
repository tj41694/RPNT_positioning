# CMakeList.txt: GNSSNav 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-std=c++11)
find_package(Threads)

file(GLOB SRC_FILES
  "postmain.cpp"
  "./BaseFunction/*.cpp"
  "./Adjustment/*.cpp"
  "./Decode/*.cpp"
  "./Decode/raw/*.cpp"
  "./GNSS/*.cpp"
  "./GNSS/AmbModel/*.cpp"
  "./GNSS/AntModel/*.cpp"
  "./GNSS/DataClass/*.cpp"
  "./GNSS/EphModel/*.cpp"
  "./GNSS/IonModel/*.cpp"
  "./GNSS/ParModel/*.cpp"
  "./GNSS/PosModel/*.cpp"
  "./GNSS/TidModel/*.cpp"
  "./GNSS/TroModel/*.cpp"
  "./InputFile/*.cpp"
  "./RtkStream/*.cpp"
  "./ConfigFile/*.cpp")

include_directories ("./")

# 将源代码添加到此项目的可执行文件。
add_executable (GNSSNav ${SRC_FILES})
target_link_libraries(GNSSNav ${CMAKE_THREAD_LIBS_INIT})

# TODO: 如有需要，请添加测试并安装目标。
