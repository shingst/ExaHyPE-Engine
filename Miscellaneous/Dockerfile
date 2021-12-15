FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y python3 libtbb-dev git build-essential cmake mpich libnetcdf-dev libnuma-dev doxygen libyaml-cpp-dev libeigen3-dev

RUN cd / && mkdir dependencies && \
    cd /dependencies/ && git clone https://github.com/uphoffc/ImpalaJIT.git && \
    cd ImpalaJIT/ && \
    mkdir build && cd build &&\
    cmake -DSHARED_LIB=1 .. && make && make install

ENV PKG_CONFIG_PATH="/dependencies/ImpalaJIT/build/:/usr/lib/x86_64-linux-gnu/pkgconfig"

RUN cd /dependencies && git clone --recursive https://github.com/TUM-I5/ASAGI.git && \
    cd ASAGI/ && \
    mkdir build && cd build &&\
    cmake .. && make && make install

ENV COMPILER_LFLAGS=" -L/dependencies/ImpalaJIT/build/ -limpalajit -L/usr/lib/x86_64-linux-gnu/ -lyaml-cpp  -L/dependencies/ASAGI/ -lasagi -lpthread"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/dependencies/ASAGI/build:/dependencies/ImpalaJIT/build"

RUN cd /dependencies/ && git clone https://github.com/SeisSol/easi.git && \
    cd easi && git checkout 18382bf60204c67782057fc371c1e699c9bb31b0 && \
    sed -i "s|find_package (YAMLCPP REQUIRED yaml-cpp)|pkg_check_modules (YAMLCPP REQUIRED yaml-cpp) |g" CMakeLists.txt &&\
    sed -i "s|find_package (OpenMP)| |g" CMakeLists.txt && \
    mkdir build && cd build && CC=mpicc CXX=mpicxx cmake .. && make

RUN git clone --recurse-submodules https://gitlab.lrz.de/exahype/ExaHyPE-Engine.git
WORKDIR /ExaHyPE-Engine
RUN git checkout ci_dev && \ 
    git pull
RUN Submodules/updateSubmodules.sh

ENV COMPILER_LFLAGS="  -L/dependencies/easi/build/ -leasi -L/dependencies/ImpalaJIT/build/ -limpalajit -L/usr/lib/x86_64-linux-gnu/ -lyaml-cpp  -L/dependencies/ASAGI/ -lasagi -lpthread"
ENV COMPILER_CFLAGS=" -I/dependencies/easi/include/"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/dependencies/easi/build"
ENV COMPILER=GNU
ENV TBB_INC=/usr/include/tbb
ENV TBB_SHLIB=/usr/lib/x86_64-linux-gnu/libtbb.so
ENV EXAHYPE_CC=mpicc
ENV EXAHYPE_CXX=mpicxx
ENV PROJECT_CFLAGS=-DUSE_ASAGI

#RUN cd ApplicationExamples/TestCases/nonlinear_Euler_Sod_3d && \
#    ../../../Toolkit/toolkit.sh nonlinear_Euler_Sod_3d.exahype2 && \
#    make -j2

RUN cd ApplicationExamples/TestCases/nonlinear_Euler_Sod_3d && \
    ../../../Toolkit/toolkit.sh nonlinear_Euler_Sod_3d_generic.exahype2 && \
    make -j2
