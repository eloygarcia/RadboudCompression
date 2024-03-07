# FROM nvidia/cuda:11.7.1-devel-ubuntu20.04
FROM ubuntu:20.04
LABEL authors="eloygarcia"

#ENTRYPOINT ["top", "-b"]

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
# RUN apt-get install -y tzdata

RUN set -ex;		\
apt-get update; 	\
apt-get install -y g++ curl libzmq3-dev wget git; \
apt-get install -y libx11-dev libxt-dev xserver-xorg-dev xorg-dev ;

# RUN apt-get update; apt-get install -y build-essential
# Install OpenGL
# Debian, Ubuntu
# https://en.wikibooks.org/wiki/OpenGL_Programming/Installation/Linux
RUN apt-get update && apt-get install --yes build-essential libgl1-mesa-dev libssl-dev

## Intalling cgal
RUN  apt-get install -y libcgal-dev

# cmake
RUN wget http://www.cmake.org/files/v3.16/cmake-3.16.6.tar.gz \
    && tar xzf cmake-3.16.6.tar.gz
# old version 2.8.12.2
# RUN wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz \
#    && tar xzf cmake-2.8.12.2.tar.gz

RUN cd cmake-3.16.6 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install


# cmake
#RUN apt-get update && apt-get -y install cmake protobuf-compiler cmake-gui

# Download & build Tcls
# https://www.tcl.tk/doc/howto/compile.html#unix
RUN wget http://prdownloads.sourceforge.net/tcl/tcl8.6.6-src.tar.gz && tar -zxvf tcl8.6.6-src.tar.gz
RUN cd tcl8.6.6/unix && ./configure && make && make install
#
# Download & build Tk
# https://www.tcl.tk/doc/howto/compile.html
RUN wget http://prdownloads.sourceforge.net/tcl/tk8.6.6-src.tar.gz && tar -zxvf tk8.6.6-src.tar.gz
RUN cd tk8.6.6/unix && ./configure && make && make install

# vtk
# Download & extract VTK
# check : https://discourse.vtk.org/t/installing-vtk-in-ubuntu-18-04/2147/3
RUN wget https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz && tar -zxvf VTK-9.3.0.tar.gz

# Required version 9.3.0 or higher due to "vtkStaticCleanUnstructuredGrid.h"


RUN mkdir -p /home/vtk-build
RUN cd /home/vtk-build/ && cmake \
   -DCMAKE_BUILD_TYPE:STRING=Release \
   -DBUILD_TESTING:BOOL=OFF \
   -DVTK_WRAP_PYTHON:BOOL=OFF \
   -DVTK_WRAP_PYTHON_SIP:BOOL=OFF \
   -DVTK_WRAP_TCL:BOOL=ON \
   -DVTK_PYTHON_VERSION:STRING=3 \
   -DVTK_USE_TK:BOOL=ON \
   /VTK-9.3.0/ && make

## Download & build itk \
RUN wget https://github.com/InsightSoftwareConsortium/ITK/releases/download/v5.3.0/InsightToolkit-5.3.0.tar.gz && tar -zxvf InsightToolkit-5.3.0.tar.gz

RUN mkdir -p /home/itk-build
RUN cd /home/itk-build && cmake \
    -DBUILD_EXAMPLES:BOOL=OFF \
    -DITK_BUILD_DEFAULT_MODULES:BOOL=ON \
    -DModule_IOMeshSWC:BOOL=ON \
    -DModule_ITKVtkGlue:BOOL=ON \
    -DVTK_DIR:PATH=/home/vtk-build \
    -DZLIBNG_ENABLE_TESTS:BOOL=ON \
    /InsightToolkit-5.3.0/ && make


## Copying source code
COPY ./ $HOME/home/Model-Extraction-source
RUN mkdir -p /home/Model-Extraction
RUN cd /home/Model-Extraction && cmake \
    -G "Unix Makefiles" \
    -DCMAKE_BUILD_TYPE:STRING=Release -DCGAL_ImageIO_USE_ZLIB:BOOL=ON \
    -DCGAL_DIR:STRING=/usr/lib/x86_64-linux-gnu/cmake/CGAL \
    -DITK_DIR:STRING=/home/itk-build \
    -DVTK_DIR:STRING=/home/vtk-build \
    /home/Model-Extraction-source/ && make
