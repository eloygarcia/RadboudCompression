cmake_minimum_required(VERSION 3.5)
project(ItkToCGAL)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -frounding-math")
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

# find_package(ItkVtkGlue REQUIRED)
# include(${ItkVtkGlue_USE_FILE})
# set(Glue ItkVtkGlue)
 
find_package(CGAL REQUIRED COMPONENTS Core ImageIO) # GMP RS3 MPFI )
# include_directories(${CGAL_INCLUDE_DIRS})
# set(LIBS ${LIBS} ${CGAL_LIBRARIES})
include(${CGAL_USE_FILE})

find_package( Boost REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS
    ${LIBS}
    ${Boost_LIBRARIES})

set(SOURCE_FILES source.cpp)
add_executable(${PROJECT_NAME}
                ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME}
                         ${CGAL_LIBS}
                         ${Boost_LIBRARIES}
                         ${Glue}
                         ${VTK_LIBRARIES}
                         ${ITK_LIBRARIES})
