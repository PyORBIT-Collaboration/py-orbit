[![Build Status](https://travis-ci.org/PyORBIT-Collaboration/py-orbit.svg?branch=master)](https://travis-ci.org/PyORBIT-Collaboration/py-orbit)

# Py-ORBIT  Installation

Installation procedure requires building from source.
All installation steps happen in command line (terminal).

## 1. Installing required libraries

### Ubuntu (and other distributions using apt-get: Debian, Mint etc)
```shell
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install python-dev libmpich-dev mpich  zlib1g-dev libfftw3-dev
```

### RedHat (and other distributions using yum: Fedora, CentOS etc)
```shell
sudo yum update
sudo yum group install "Development Tools"
sudo yum install python-devel mpich mpich-devel zlib-devel fftw-devel
```

### Mac 
There are different package managers available.

[MacPorts](https://www.macports.org)  
Sync with package repository

`sudo port -v selfupdate`

After syncing run:
```shell
sudo port install fftw mpich
```
Alternatively you can use [Homebrew](http://brew.sh). 
This tends to build a lot of packages from sources (especially for the first time), which can take a long time.
```shell
brew update
brew install fftw mpich
```

### Building the whole environment from source
If you don't want to use standard libraries supplied by your distribution, you can build the whole environment from scratch. It is also possible to do this without having a root account. The process is described in detail [here](BuildFromSource.md).

## 2. Clone the source code
```shell
git clone https://github.com/PyORBIT-Collaboration/py-orbit.git
```
Your source is now in the *py-orbit* directory.
## 3. Setup environment variables
*setupEnevironment.sh* will try to figure out all paths.
This should be sufficient for common Linux distributions. If you built the environment form source, use *customEnvironment.sh* instead.
```shell
cd py-orbit
source setupEnvironment.sh
```


## 4. Build the code

```shell 
make clean
make
```
If make failed, it usually means that some of the libraries aren't set up properly.



# Running Examples

Setup the environment variables (needs to be done once per teminal session). If you built the environment from source, use *customEnvironment.sh* instead.
Alternatively you can place `source <path-to-pyORBIT-installation>/setupEnvironment.sh` in your *.bashrc*.
```shell
source setupEnvironment.sh
cd examples/AccLattice_Tests
./START.sh lattice_test.py 2
```
This will launch *lattice_test* example on two MPI nodes. Other examples are availabale in [Examples](../../../examples) repository.

# Structure
**./src**		- source code for the core ORBIT C++ classes, including
		  wrappers, etc.

**./py**		- python modules and wrapper classes for the core ORBIT
		  classes.

**./ext**		- source code for external modules. Compilations of this
		  code should be placed into **./lib**.

**./lib**  	- .so shared libraries to be used under pyORBIT interpreter.

**./doc**		- pyORBIT documentation.

**./examples**		- pyORBIT examples.

**./conf**		- configuration information.

**./bin**		-  pyORBIT executables.
