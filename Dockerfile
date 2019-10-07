FROM centos
WORKDIR /pyORBIT
ADD . /pyORBIT

# RUN yum update
RUN yum install -y gcc gcc-c++ make python2 python2-devel mpich mpich-devel zlib-devel fftw-devel
RUN bash -c "cd /pyORBIT; source setupEnvironment.sh; make clean; make"
