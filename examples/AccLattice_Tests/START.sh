#!/bin/bash

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the python script> <N-CPUs>"
    exit $E_BADARGS
fi

if [ ! -n "$2" ]
  then
    echo "Usage: `basename $0` <name of the python script> <N CPUs>"
    exit $E_BADARGS
fi

mpirun -np $2 ${ORBIT_ROOT}/bin/pyORBIT $1
