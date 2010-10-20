import math

from bunch import Bunch

def bunch_orbit_to_pyorbit(ringLength, kineticEnergy, name_of_orbit_mpi_bunch_file, pyOrbitBunch = None):
	"""
	Translates ORBIT_MPI bunch to pyORBIT bunch and returns it. PyORBIT bunch needs 
	the ring length (m) and energy, mass and charge of the synchronous particle, but 
	the ORBIT_MPI does not have it. So, this information is specified in pyOrbitBunch or
	it will be proton by default.
	ORBIT_MPI file has lines: x[mm] xp[mrad] y[mm] yp[mrad]   phi[rad]  dE[GeV].
	pyORBIT: x[m] xp[rad] y[m] yp[rad]  z[m]  dE[GeV]
	"""
	L =  ringLength
	file_in = open(name_of_orbit_mpi_bunch_file,"r")
	if(pyOrbitBunch == None):  pyOrbitBunch = Bunch()
	pyOrbitBunch.getSyncParticle().kinEnergy(kineticEnergy)

	ln = file_in.readline()
	while ln:
		res_arr = ln.strip().split()
		val_arr = []
		for s in res_arr:
			val_arr.append(float(s))
		for i in range(4):
			val_arr[i] /= 1000.
		val_arr[4] = val_arr[4]*L/(2*math.pi)
		val_arr[5] = val_arr[5]
		pyOrbitBunch.addParticle(val_arr[0],val_arr[1],val_arr[2],val_arr[3],val_arr[4],val_arr[5])
		ln = file_in.readline()
		
	file_in.close()
	return pyOrbitBunch

