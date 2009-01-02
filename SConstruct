import os
from distutils import sysconfig

path = os.environ["PATH"]

mpi_cpp = os.environ["MPI_CPP"]

cpp_flags = ""
if os.environ.has_key("CCFLAGS"):
	cpp_flags = os.environ["CCFLAGS"]

cpp_flags = cpp_flags + " -fPIC " 
cpp_flags = cpp_flags + " -DUSE_MPI=1 -DMPICH_IGNORE_CXX_SEEK " 

cpp_shared_lib_flags = " -Xlinker -export-dynamic "

Decider('make')

#--------------------------------------
# Make pyORBIT
#--------------------------------------
incl_dirs = []
incl_dirs.append("./src/main/")
incl_dirs.append("./src/mpi/")
incl_dirs.append("./src/orbit/")
incl_dirs.append("./src/teapot/")
incl_dirs.append("./src/trackerrk4/")
incl_dirs.append("./src/utils/")

py_lib_path = sysconfig.get_config_var('LIBPL')
py_libs = [sysconfig.get_config_var('LIBRARY')]
py_shared_libs = sysconfig.get_config_var('SHLIBS').split()
py_include_dir = sysconfig.get_config_var('INCLUDEPY')

pyOrbit_incl_dirs = incl_dirs + [py_include_dir,]
py_libs = py_libs + py_shared_libs

pyOrbitEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, LINKFLAGS =cpp_shared_lib_flags,  CPPPATH = pyOrbit_incl_dirs,  ENV = {'PATH':path})

cpp_files_list = []
for dr in incl_dirs:
	pyOrbitEnv.VariantDir(dr+"obj", dr, duplicate=0)
	cpp_files_list = cpp_files_list + Glob(dr+"obj"+"/*.cc")

pyOrbit_exe = pyOrbitEnv.Program('./bin/pyORBIT',cpp_files_list, LIBS = py_libs , LIBPATH = py_lib_path)
Default(pyOrbit_exe)


#--------------------------------------
# Make extensions 
#--------------------------------------

#--------- Laser Stripping Module -----
incl_ext_dirs = []
incl_ext_dirs.append("./ext/laserstripping/")

incl_dirs = pyOrbit_incl_dirs + incl_ext_dirs
laserStrippingFieldEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, CPPPATH = incl_dirs,  ENV = {'PATH':path})

cpp_files_list = []
for dr in incl_ext_dirs:
	laserStrippingFieldEnv.VariantDir(dr+"obj", dr, duplicate=0)
	cpp_files_list = cpp_files_list + Glob(dr+"obj"+"/*.cc")

laserstripping_lib = laserStrippingFieldEnv.SharedLibrary('./lib/laserstripping',
	                          cpp_files_list, 
														#LIBS = py_libs,
                            #LIBPATH = py_lib_path, 
                            LINKFLAGS = cpp_shared_lib_flags,
														SHLIBPREFIX = "")
Default(laserstripping_lib)

#--------- Space Charge 2D Field Module -----
incl_ext_dirs = []
incl_ext_dirs.append("./ext/spacecharge/")
if (not os.environ.has_key("FFTW3_ROOT")):
	print "You have to define FFTW3_ROOT env variable"
	sys.exit(1)
incl_ext_dirs.append(os.environ["FFTW3_ROOT"]+"/include")

incl_dirs = pyOrbit_incl_dirs + incl_ext_dirs
spacechargeEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, CPPPATH = incl_dirs,  ENV = {'PATH':path})

cpp_files_list = []
for dr in incl_ext_dirs:
	spacechargeEnv.VariantDir(dr+"obj", dr, duplicate=0)
	cpp_files_list = cpp_files_list + Glob(dr+"obj"+"/*.cc")

spacecharge_lib = spacechargeEnv.SharedLibrary('./lib/spacecharge',
	                          cpp_files_list, 
														LIBS = ["libfftw3"],
                            LIBPATH = [os.environ["FFTW3_ROOT"]+"/lib",] ,
                            LINKFLAGS = cpp_shared_lib_flags,
														SHLIBPREFIX = "")
Default(spacecharge_lib)

#--------------------------------------
# Make documentation (see Phony Target for Scons)
#--------------------------------------
docEnv = Environment()
docEnv.AlwaysBuild(docEnv.Alias('docs',[],"doxygen"))

