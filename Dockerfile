FROM centos
WORKDIR /pyORBIT
ADD . /pyORBIT

# RUN yum update
RUN yum install -y gcc gcc-c++ make python-devel mpich mpich-devel zlib-devel fftw-devel
RUN bash -c "cd pyORBIT; source setupEnvironment.sh; make clean; make"
