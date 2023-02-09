cmake_minimum_required(VERSION 3.5)
project(WriteNiftyModel)

find_package(ITK)
include(${ITK_USE_FILE})

find_package(VTK)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

add_library( Library xml.cpp
                     xmlModelWriter.cpp
					 NiftySimEjecutable.cpp)

add_executable( WriteNiftyModel Source.cpp)

target_link_libraries(	WriteNiftyModel
						${VTK_LIBRARIES}
						${ITK_LIBRARIES}
						Library )