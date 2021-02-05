FROM ubuntu:bionic
CMD ["/bin/sh"]
MAINTAINER FLB

ARG TMPDIR=/tmp
RUN apt-get update \
 && export DEBIAN_FRONTEND=noninteractive \
 && apt-get -y install bash vim sed wget sudo git\
     make cmake g++ gcc gfortran mpich \
     libmpich-dev libfftw3-dev libblas-dev liblapack-dev libhdf5-dev libnetcdf-c++4-dev \ 
     python3 python3-numpy python3-netcdf4 python3-scipy python3-matplotlib ipython3 python3-pip\
 && sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.6 1 \
 && sudo update-alternatives --install /usr/bin/ipython ipython /usr/bin/ipython3 1 \
 && apt-get clean

# Install sundials
WORKDIR /
RUN wget https://computation.llnl.gov/projects/sundials/download/sundials-4.1.0.tar.gz \
 && tar -xvzf sundials-4.1.0.tar.gz \
 && cd sundials-4.1.0  \
 && mkdir build \
 && mkdir -p /usr/examples \
 && cd build \        
 && cmake -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -DBUILD_STATIC_LIBS:BOOL=true \
           -DBUILD_SHARED_LIBS:BOOL=true -DEXAMPLES_ENABLE:BOOL=true \ 
           -DEXAMPLES_INSTALL:BOOL=true -DMPI_ENABLE:BOOL=true \
           -DOPENMP_ENABLE:BOOL=true -DPTHREAD_ENABLE:BOOL=true \
           -DCUDA_ENABLE:BOOL=false -DRAJA_ENABLE:BOOL=false -DLAPACK_ENABLE=ON \
           -DFCMIX_ENABLE=ON -DCMAKE_C_FLAGS="-fPIC" -DSUNDIALS_INDEX_SIZE=32 ..\  
 && make \
 && make install \
 && cd / \
 && rm -rf sundials-4.1.0 \

RUN cd /home
RUN pip3 install sympy ase
#RUN git clone https://github.com/francescalb/micki.git 
#RUN export PYTHONPATH=$HOME/micki:$PYTHONPATH
WORKDIR /home
