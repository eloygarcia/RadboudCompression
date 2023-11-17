cmake_minimum_required(VERSION 3.5)
project(WriteNiftyModel)

find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

# find_package(ITK REQUIRED)
# include(${ITK_USE_FILE})

#add_library( Library x)

add_executable( WriteNiftyModel Source.cpp
                     xml.cpp
                     xmlModelWriter.cpp
					 NiftySimEjecutable.cpp)

target_link_libraries(	WriteNiftyModel
						${VTK_LIBRARIES} )
						# ${ITK_LIBRARIES} )
						# Library )