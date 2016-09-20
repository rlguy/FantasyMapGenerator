cmake_minimum_required(VERSION 2.8.9)

project(map_generation)
set(CMAKE_BUILD_TYPE Release)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W4")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall")
endif()

file(GLOB SOURCES "src/*.cpp")

add_library(objects OBJECT ${SOURCES})

add_executable(map_generation $<TARGET_OBJECTS:objects>)

file(COPY "src/font_data" DESTINATION ${CMAKE_BINARY_DIR})
file(COPY "src/city_data" DESTINATION ${CMAKE_BINARY_DIR})