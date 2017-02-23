cmake_minimum_required(VERSION 2.8.9)

project(map_generation)
set(CMAKE_BUILD_TYPE Release)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W4 /D_CRT_SECURE_NO_WARNINGS")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -Doff64_t=_off64_t")
endif()

find_package(PythonLibs)
if (PYTHONLIBS_FOUND)

	if(${PYTHONLIBS_VERSION_STRING} VERSION_LESS 3)
		add_definitions(-DPYTHON_VERSION_2)
	else()
		add_definitions(-DPYTHON_VERSION_3)
	endif()

	set(CHECK_CAIRO_PATH "${CMAKE_SOURCE_DIR}/src/render/check_cairo_support.py")
	execute_process(COMMAND python ${CHECK_CAIRO_PATH} RESULT_VARIABLE IS_CAIRO_SUPPORTED)
	if (IS_CAIRO_SUPPORTED)
		add_definitions(-DPYTHON_RENDERING_SUPPORTED)
	else()
		message(WARNING "Warning: Pycairo was not detected on target machine. The Python cairo graphics library (http://cairographics.org/pycairo/) is required for map drawing.\nProject will be built without drawing capabilities.")
	endif()
else()
	message(WARNING "Warning: Python was not detected on target machine. Python and the Pycairo graphics library are required for map drawing.\nProject will be built without drawing capabilities.")
endif()

include_directories(src ${PYTHON_INCLUDE_DIRS})
file(GLOB SOURCES "src/*.cpp" "src/argtable3/*.c")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})

add_library(objects OBJECT ${SOURCES})
add_executable(map_generation $<TARGET_OBJECTS:objects>)
target_link_libraries(map_generation ${PYTHON_LIBRARIES})

file(COPY "src/fontdata" DESTINATION ${CMAKE_BINARY_DIR})
file(COPY "src/citydata" DESTINATION ${CMAKE_BINARY_DIR})
file(COPY "src/render" DESTINATION "${CMAKE_BINARY_DIR}")

set(RESOURCES_EXECUTABLE_DIRECTORY  ${CMAKE_BINARY_DIR})
set(RESOURCES_FONT_DATA_DIRECTORY   ${CMAKE_BINARY_DIR}/fontdata)
set(RESOURCES_CITY_DATA_DIRECTORY   ${CMAKE_BINARY_DIR}/citydata)
set(RESOURCES_FONT_DATA_RESOURCE    ${CMAKE_BINARY_DIR}/fontdata/fontdata.json)
set(RESOURCES_CITY_DATA_RESOURCE    ${CMAKE_BINARY_DIR}/citydata/countrycities.json)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/resources.h.in"
  "${PROJECT_SOURCE_DIR}/src/resources.h"
)