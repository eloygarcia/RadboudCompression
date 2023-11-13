cmake_minimum_required(VERSION 3.5)
project(FeBioModel)

find_package(VTK)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS  "-std=c++0x ${CMAKE_CXX_FLAGS}")

add_executable( FeBioModel paddleConstruction.cpp )

find_package(VTK COMPONENTS
  CommonColor
  CommonCore
  CommonDataModel
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)
target_link_libraries(	FeBioModel
						${VTK_LIBRARIES} )