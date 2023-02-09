cmake_minimum_required(VERSION 3.5)
project(ReconstructImage)

find_package(ITK)
include(${ITK_USE_FILE})

find_package(VTK)
# include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

# add_library( Library RegularGrid.cpp NewProjection.cpp)

find_package(CUDA REQUIRED)
# include_directories(
# 						${CUDA_INCLUDE_DIRS}
# 						${CUDA_SDK_COMMON_INCLUDE_DIR}
# 						)
# 
cuda_include_directories(	
							${CUDA_INCLUDE_DIRS}	
							${CUDA_SDK_COMMON_INCLUDE_DIR}
							)

# find_path(CUDA_SDK_COMMON_INCLUDE_DIR
#       helper_cuda.h
#       PATHS ${CUDA_SDK_SEARCH_PATH}
#      PATH_SUFFIXES "common/inc"
#       DOC "Location of helper_cuda.h"
#     NO_DEFAULT_PATH
 #      )

# CUDA_INCLUDE_DIRECTORIES(
# 							${CMAKE_SOURCE_DIR}/src/
# 							${CMAKE_SOURCE_DIR}/source
#							)
# cuda_compile(HELPER_O helper.cu)
# cuda_compile(DRIVER_O driver.cu OPTIONS -G)
 
# set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ;-D_GPU_" )
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GPU_" )
#  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --gpu-architecture sm_20)
# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} --gpu-architecture sm_20 )
# set(CUDA_SEPARABLE_COMPILATION OFF)  # En off evita un error de compilacion del cuda !

#set(CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" )
#include(FindCUDA)

cuda_add_library( CudaLibrary RegularGrid.cpp
							CudaProjection.cu
							CudaLibrary_kernels.cu
							NewProjection.cpp )

add_executable( ReconstructImage Source.cpp)

target_link_libraries(	ReconstructImage
						${VTK_LIBRARIES}
                        ${ITK_LIBRARIES} 
						CudaLibrary )