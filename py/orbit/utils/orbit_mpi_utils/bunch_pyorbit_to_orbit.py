import math

from bunch import Bunch

def bunch_pyorbit_to_orbit(ringLength, pyOrbitBunch, name_of_orbit_mpi_bunch_file):
	"""
	Translates pyORBIT bunch to ORBIT_MPI bunch and dumps it into the file.
	The ring length should be defined in the input (in meters).
	ORBIT_MPI file has lines: x[mm] xp[mrad] y[mm] yp[mrad]   phi[rad]  dE[GeV].
	pyORBIT: x[m] xp[rad] y[m] yp[rad]  z[m]  dE[GeV]
	"""	
	L = ringLength
	file_out = open(name_of_orbit_mpi_bunch_file,"w")
	b = pyOrbitBunch
	for i in range(b.getSize()):
		x = b.x(i)*1000.
		px = b.px(i)*1000.
		y = b.y(i)*1000.
		py = b.py(i)*1000.
		z = b.z(i)*2.0*math.pi/L
		dE = b.dE(i)
		file_out.write(str(x) + " " + str(px) + " " + str(y) + " " + str(py) + " "+ str(z) + " " + str(dE) + "\n")	
	file_out.close()

