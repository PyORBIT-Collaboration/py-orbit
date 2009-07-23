import os
from distutils import sysconfig

path = os.environ["PATH"]

mpi_cpp = os.environ["MPI_CPP"]
	
cpp_flags = ""
if os.environ.has_key("CXXFLAGS"):
	cpp_flags = os.environ["CXXFLAGS"]
	
cpp_link_flags =""
if os.environ.has_key("LINKFLAGS"):
	 cpp_link_flags = os.environ["LINKFLAGS"]
	
cpp_flags = cpp_flags + " -fPIC " 
cpp_flags = cpp_flags + " -DUSE_MPI=1 -DMPICH_IGNORE_CXX_SEEK " 

cpp_shared_lib_flags = " -Xlinker -export-dynamic "

Decider('make')

#-----We need FFTW3 library
if (not os.environ.has_key("FFTW3_ROOT")):
	print "You have to define FFTW3_ROOT env variable"
	sys.exit(1)

#---Function returns the list of 	sub_directories of d_path
def getSubDirs(d_path,ignore_list):
	dirs = [f for f in os.listdir(d_path) if os.path.isdir(os.path.join(d_path,f))]
	return filter(lambda x: not x in ignore_list, dirs)

#---Function returns the dictionary {key=sub_dir:val=[sub_sub_dirs,...]} 
#---   sub_dir is a	subdirectory of d_path
def getDictSubDirs(d_path,ignore_list):
	s_dirs = getSubDirs(d_path,ignore_list)
	s_dirs_dict = {}
	for dr in s_dirs:
		s_dirs_dict[dr] = getSubDirs(os.path.join(d_path,dr),ignore_list)
	return s_dirs_dict
	
#---Function returns the list with INCLUDE directories for compilation
def getInclSubDirs(d_path,ignore_list):
	s_dirs = getSubDirs(d_path,ignore_list)
	s_dirs_dict = getDictSubDirs(d_path,ignore_list)
	inc_dirs = []
	for dr in s_dirs:
		inc_dirs.append(os.path.join(d_path,dr))
		for s_dr in s_dirs_dict[dr]:
			inc_dirs.append(os.path.join(d_path,dr,s_dr))
	return inc_dirs
		
#---Function will make Variants "obj" directories - env is the Environment() instance 
#--- This function returns the list of source file "*.cc" for compilation
def makeVariantDirs(env,d_path,ignore_list):
	s_dirs = getSubDirs(d_path,ignore_list)
	s_dirs_dict = getDictSubDirs(d_path,ignore_list)
	cpp_files_list = []
	for dr in s_dirs:
		varDir = os.path.join(d_path,dr,"obj")
		srcDir = os.path.join(d_path,dr)
		pyOrbitEnv.VariantDir(varDir,srcDir,duplicate=0)
		cpp_files_list = cpp_files_list + Glob(os.path.join(varDir,"*.cc"))
		for s_dr in s_dirs_dict[dr]:
			cpp_files_list = cpp_files_list + Glob(os.path.join(varDir,s_dr,"*.cc"))
	return cpp_files_list

#--------------------------------------
# Make pyORBIT
#--------------------------------------
m_path = "src"
orbit_src_dirs = getSubDirs(m_path,[".svn","tests"])
orbit_incl_dirs = getInclSubDirs(m_path,[".svn","tests","obj"])
		
orbit_incl_dirs.append(os.environ["FFTW3_ROOT"]+"/include")

py_lib_path = [sysconfig.get_config_var('LIBPL'),os.environ["FFTW3_ROOT"]+"/lib"]
py_libs = [sysconfig.get_config_var('LIBRARY')]
py_shared_libs = sysconfig.get_config_var('SHLIBS').split()
py_include_dir = sysconfig.get_config_var('INCLUDEPY')

pyOrbit_incl_dirs = orbit_incl_dirs + [py_include_dir,]
py_libs = py_libs + py_shared_libs
py_libs = py_libs + ["libfftw3",]


pyOrbitEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, \
						 LINKFLAGS =cpp_shared_lib_flags + cpp_link_flags,\
						 CPPPATH = pyOrbit_incl_dirs,  ENV = {'PATH':path})

cpp_files_list = makeVariantDirs(pyOrbitEnv,m_path,[".svn","tests","obj"])

pyOrbit_exe = pyOrbitEnv.Program('./bin/pyORBIT',cpp_files_list, LIBS = py_libs , LIBPATH = py_lib_path)

Default(pyOrbit_exe)

#--------------------------------------
# Make extensions 
#--------------------------------------
ext_dir = "ext"

#--------- Laser Stripping Module -----
project_dir = "laserstripping"
exclude_dirs = getSubDirs(ext_dir,[])
exclude_dirs.remove(project_dir)

incl_ext_dirs = getInclSubDirs(ext_dir,[".svn","obj","transitions","PyModules"]+exclude_dirs)

incl_dirs = pyOrbit_incl_dirs + incl_ext_dirs
laserStrippingFieldEnv = Environment(CXX = mpi_cpp, CCFLAGS = cpp_flags, CPPPATH = incl_dirs,  ENV = {'PATH':path})

cpp_files_list = makeVariantDirs(laserStrippingFieldEnv,ext_dir,[".svn","obj","transition","PyModules"]+exclude_dirs)

laserstripping_lib = laserStrippingFieldEnv.SharedLibrary('./lib/laserstripping',
	                          cpp_files_list, 
														#LIBS = py_libs,
                            #LIBPATH = py_lib_path, 
                            LINKFLAGS = cpp_shared_lib_flags + cpp_link_flags,
														SHLIBPREFIX = "")
Default(laserstripping_lib)

#--------------------------------------
# Make documentation (see Phony Target for Scons)
#--------------------------------------
docEnv = Environment()
docEnv.AlwaysBuild(docEnv.Alias('docs',[],"doxygen"))

