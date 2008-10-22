import os
from distutils import sysconfig

path = os.environ["PATH"]

mpi_cpp = os.environ["MPI_CPP"]

cpp_flags = ""
if os.environ.has_key("CCFLAGS"):
	cpp_flags = os.environ["CCFLAGS"]

cpp_flags = cpp_flags + " -fPIC " 
cpp_flags = cpp_flags + " -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK " 

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
incl_dirs.append("./src/utils/")

py_lib_path = sysconfig.get_config_var('LIBPL')
py_libs = [sysconfig.get_config_var('LIBRARY')]
py_shared_libs = sysconfig.get_config_var('SHLIBS').split()
py_include_dir = sysconfig.get_config_var('INCLUDEPY')

pyOrbit_incl_dirs = incl_dirs + [py_include_dir,]
py_libs = py_libs + py_shared_libs

pyOrbitEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, CPPPATH = pyOrbit_incl_dirs,  ENV = {'PATH':path})

cpp_files_list = []
for dr in incl_dirs:
	pyOrbitEnv.VariantDir(dr+"obj", dr, duplicate=0)
	cpp_files_list = cpp_files_list + Glob(dr+"obj"+"/*.cc")

pyOrbit_exe = pyOrbitEnv.Program('./bin/pyORBIT',cpp_files_list, LIBS = py_libs , LIBPATH = py_lib_path)
Default(pyOrbit_exe)


#--------------------------------------
# Make extensions 
#--------------------------------------

#--------- Traker 3D Field Module -----
incl_ext_dirs = []
incl_ext_dirs.append("./ext/tracker3dfield/")

incl_dirs = pyOrbit_incl_dirs + incl_ext_dirs
tracker3DFieldEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, CPPPATH = incl_dirs,  ENV = {'PATH':path})

cpp_files_list = []
for dr in incl_ext_dirs:
	tracker3DFieldEnv.VariantDir(dr+"obj", dr, duplicate=0)
	cpp_files_list = cpp_files_list + Glob(dr+"obj"+"/*.cc")

tracker3D_lib = tracker3DFieldEnv.SharedLibrary('./modules/tracker3dfield',
	                          cpp_files_list, LIBS = py_libs,
                                  LIBPATH = py_lib_path, 
                                  LINKFLAGS = cpp_shared_lib_flags)
Default(tracker3D_lib)

#--------------------------------------
# Make documentation
#--------------------------------------
if 'docs' in COMMAND_LINE_TARGETS:
	import posix
	posix.system("doxygen")

#---------------------------------------
# Clean documentation
#---------------------------------------
print "COMMAND_LINE_TARGETS =",COMMAND_LINE_TARGETS
if '-c' in COMMAND_LINE_TARGETS:
	print "DEBUG"	
