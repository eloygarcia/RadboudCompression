cmake_minimum_required(VERSION 2.8)
project(WriteNiftyModel)

find_package(VTK)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

add_library( Library xml.cpp
                     xmlModelWriter.cpp)

add_executable( WriteNiftyModel Source.cpp)

target_link_libraries(	WriteNiftyModel
						${VTK_LIBRARIES} 
						Library )