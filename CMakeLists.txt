cmake_minimum_required(VERSION 3.23)
project(Demo)

set(CMAKE_CXX_STANDARD 14)

include_directories(MyTF)

add_executable(Demo main.cpp)
