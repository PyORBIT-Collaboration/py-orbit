# Building required libraries.
All libraries are built from source codes and installed in specific locations.
You need to download sources first. The following versions are known to work, but you can always try out a different version.
```shell
mkdir pyORBIT
cd pyORBIT
```
Everything will be contained in the pyORBIT directory.

## Download all required sources and extract them.
```shell
curl http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz | tar xvz
curl https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz | tar xvz
curl http://zlib.net/fossils/zlib-1.2.11.tar.gz | tar xvz
curl http://www.fftw.org/fftw-3.3.5.tar.gz | tar xvz
```

## Build all libraries
All libraries will be installed in the local directory (pyORBIT).

1. Build Python
```shell
cd Python-2.7.12
./configure -prefix=`pwd`/..
make
make install
cd ..
```

2. Build ZLIB
```shell
cd zlib-1.2.11
./configure -prefix=`pwd`/..
make
make install
cd ..
```

3. Build MPI
```shell
cd mpich-3.2
./configure -prefix=`pwd`/.. --disable-fortran
make
make install
cd ..
```

4. Build FFTW  (with MPI support)
```shell
cd fftw-3.3.5
./configure -prefix=`pwd`/.. --disable-fortran --enable-mpi MPICC=`pwd`/../bin/mpicc
make
make install
cd ..
```


## Proceed with installation steps
You will need to use *customEnvironment.sh* instead of *setupEnvironment.sh*.

[Step 2](README.md#2-clone-the-source-code)
