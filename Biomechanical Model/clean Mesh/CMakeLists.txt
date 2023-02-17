cmake_minimum_required(VERSION 3.5)
project(cleaningMesh)

find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

set(SOURCE_FILES Source.cpp)
add_executable(${PROJECT_NAME}
                ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME}
                         ${VTK_LIBRARIES})