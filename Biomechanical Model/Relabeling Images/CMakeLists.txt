cmake_minimum_required(VERSION 3.5)
project(RelabelingImages)

find_package(ITK)
include(${ITK_USE_FILE})

add_executable( RelabelingImages source.cpp)

target_link_libraries(	RelabelingImages
                        ${ITK_LIBRARIES} )