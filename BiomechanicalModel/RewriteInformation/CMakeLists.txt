cmake_minimum_required(VERSION 3.2)

project( RewriteInformation)

find_package(ITK)
include(${ITK_USE_FILE})

add_executable( RewriteInformation Source.cpp)
target_link_libraries( RewriteInformation ${ITK_LIBRARIES})