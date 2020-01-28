command_exists () {
    type "$1" &> /dev/null ;
}


export ORBIT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "ORBIT installed in $ORBIT_ROOT"

export ORBIT_ARCH=`uname -s`

# alias python='python2'

if command_exists python2; then
   PYEX=python2
else
   PYEX=python
fi

export PYTHON_VERSION=`$PYEX -c "from distutils import sysconfig; print sysconfig.get_config_var('VERSION');"`
echo "Python version is $PYTHON_VERSION"

PYTHON_LIB_DIR=`$PYEX -c "from distutils import sysconfig; print sysconfig.get_config_var('LIBPL');"`
if [ -f $PYTHON_LIB_DIR/libpython${PYTHON_VERSION}.a ]
   then
	export PYTHON_ROOT_LIB=$PYTHON_LIB_DIR/libpython${PYTHON_VERSION}.a
	LIB_TYPE=static
   else
	export PYTHON_ROOT_LIB="-L $PYTHON_LIB_DIR -lpython${PYTHON_VERSION}"
	LIB_TYPE=dynamic
fi

echo "Found python library: ${PYTHON_LIB_DIR} will use $LIB_TYPE library"

export PYTHON_ROOT_INC=`$PYEX -c "from distutils import sysconfig; print sysconfig.get_config_var('INCLUDEPY');"`
echo "Found Python include directory: $PYTHON_ROOT_INC"

export PYTHONPATH=${PYTHONPATH}:${ORBIT_ROOT}/py:${ORBIT_ROOT}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ORBIT_ROOT}/lib

if command_exists mpirun ; then
   echo "Found mpirun at: `which mpirun`"
   MPI_RUN_DIR=`dirname $(which mpirun)`
else
    MPI_RUN_DIR=`dirname $(find /usr 2>/dev/null| fgrep bin/mpirun | head -n1)`
    export PATH=$PATH:$MPI_RUN_DIR
    echo "Added  $MPI_RUN_DIR to PATH"
fi

export MPI_CPP=$MPI_RUN_DIR/mpicxx
echo "MPI_CPP set to $MPI_CPP"
