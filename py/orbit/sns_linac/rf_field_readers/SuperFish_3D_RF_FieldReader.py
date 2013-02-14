import sys
import math

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from spacecharge import Grid2D

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator
from orbit.utils.fitting import PolynomialFit
from orbit_utils import Polynomial

class SuperFish_3D_RF_FieldReader:
	""" 
	This class reads the SuperFish file with the 3D axial symmetric RF field. 
	It uses z and r as variables. The file include Ez, Er, and H.
	The file will work in parallel environment.
	"""
	def __init__(self):
		# self.data_arr is a flat data array with tuples [z,r,Ez,Er,E,B]
		self.data_arr = []
		self.Zmin = 0.
		self.Zmax = 0.
		self.Rmin = 0.
		self.Rmax = 0.
		self.zSteps = 0
		self.rSteps = 0
				
	def readFile(self,file_name):
		self.data_arr = []
		self.Zmin = 0.
		self.Zmax = 0.
		self.Rmin = 0.
		self.Rmax = 0.
		self.zSteps = 0
		self.rSteps = 0			
		rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
		main_rank = 0
		if(rank == 0):
			fl_in = open(file_name,"r")
			start_data = 0
			for ln in fl_in:
				res = ln.split()
				if(start_data == 0):
					if(ln.find("(Zmin,Rmin)") >= 0):
						zr_min_max = res[2][1:len(res[2])-1].split(",")
						self.Zmin = float(zr_min_max[0])
						self.Rmin = float(zr_min_max[1])
					if(ln.find("(Zmax,Rmax)") >= 0):
						zr_min_max = res[2][1:len(res[2])-1].split(",")
						self.Zmax = float(zr_min_max[0])
						self.Rmax = float(zr_min_max[1])	
					if(len(res) > 4 and res[0] == "Z" and res[1] == "and" and res[2] == "R"):
						self.zSteps = int(res[4])
						self.rSteps = int(res[5])
					if(len(res) == 6 and res[0] == '(cm)' and res[1] == '(cm)' and res[2] == '(MV/m)'):
						start_data = 1
				else:
					if(len(res) == 6):
						arr = []
						for st in res:
							arr.append(float(st))
						self.data_arr.append(arr)
					else:
						break
			fl_in.close()
		#------end of rank 0 actions
		n = len(self.data_arr)
		n = orbit_mpi.MPI_Bcast(n,mpi_datatype.MPI_INT,main_rank,mpi_comm.MPI_COMM_WORLD)
		self.zSteps = orbit_mpi.MPI_Bcast(self.zSteps,mpi_datatype.MPI_INT,main_rank,mpi_comm.MPI_COMM_WORLD)
		self.rSteps = orbit_mpi.MPI_Bcast(self.rSteps,mpi_datatype.MPI_INT,main_rank,mpi_comm.MPI_COMM_WORLD)
		self.Zmin = orbit_mpi.MPI_Bcast(self.Zmin,mpi_datatype.MPI_DOUBLE,main_rank,mpi_comm.MPI_COMM_WORLD)
		self.Zmax = orbit_mpi.MPI_Bcast(self.Zmax,mpi_datatype.MPI_DOUBLE,main_rank,mpi_comm.MPI_COMM_WORLD)
		self.Rmin = orbit_mpi.MPI_Bcast(self.Rmin,mpi_datatype.MPI_DOUBLE,main_rank,mpi_comm.MPI_COMM_WORLD)
		self.Rmax = orbit_mpi.MPI_Bcast(self.Rmax,mpi_datatype.MPI_DOUBLE,main_rank,mpi_comm.MPI_COMM_WORLD)
		if((self.zSteps+1)*(self.rSteps+1) != n):
			if(rank == 0):
				print "====================================================="
				print "SuperFish_3D_RF_FiledReader:"
				print "The file=",file_name," does not have a correct format!"
				print "Stop."
			sys.exit(1)
		for i in range(n):
			arr = [0.]*6
			if(rank == 0):
				arr = self.data_arr[i]
			arr = orbit_mpi.MPI_Bcast(arr,mpi_datatype.MPI_DOUBLE,main_rank,mpi_comm.MPI_COMM_WORLD)
			if(rank != 0):
				self.data_arr.append(arr)
		
	def getDataArray(self):
		""" 
		A convinience method. It returns the raw array 
		with records with tuples [z,r,Ez,Er,E,B]
		"""
		return self.data_arr
	
	def getNumberStepsZ(self):
		""" 
		Returns the number of steps in Z-axis. The number of grid points is self.zSteps+1.
		"""
		return self.zSteps
		
	def getNumberStepsR(self):
		""" 
		Returns the number of steps along the radius. The number of grid points is self.rSteps+1.
		"""
		return self.rSteps
	
	def makeGrid2DFileds_EzErH(self):
		""" 
		It fills out the Grid2D instances with the electric and magnetic filed components -
		Ez, Er, H.
		"""
		#The Z and R in the self.data_arr are in [cm], so to switch to [m] we use 0.0
		Zmin = 0.01*self.Zmin
		Zmax = 0.01*self.Zmax
		Rmin = 0.01*self.Rmin
		Rmax = 0.01*self.Rmax
		grid2D_Ez = Grid2D(self.zSteps+1,rSteps+1,Zmin,Zmax,Rmin,Rmax)
		grid2D_Er = Grid2D(self.zSteps+1,rSteps+1,Zmin,Zmax,Rmin,Rmax)
		grid2D_H  = Grid2D(self.zSteps+1,rSteps+1,Zmin,Zmax,Rmin,Rmax)
		for iz in range(self.zSteps+1):
			for ir in range(rSteps+1):
				i = (zSteps+1)*ir + iz
				[z,r,Ez,Er,E,B] = self.data_arr[i]
				grid2D_Ez.setValue(Ez,iz,ir)
				grid2D_Er.setValue(Er,iz,ir)
				grid2D_H.setValue(H,iz,ir)
		return (grid2D_Ez,grid2D_Er,grid2D_H)
		
	def getAxisEz(self, zSimmetric = -1):
		""" 
		Returns the Spline with Ez(z) on the axis of the RF.
		If zSimmetric > 0 the table has only half of the table,
		and the Function should be added points for (-Zmax) to (Zmin - step).
		"""
		stepZ = (self.Zmax - self.Zmin)/self.zSteps
		Ez_max = 0.
		for iz in range(self.zSteps+1):
			[z,r,Ez,Er,E,B] = self.data_arr[iz]
			Ez_abs = math.fabs(Ez)
			if(Ez_max < Ez_abs):
				Ez_max = Ez_abs			
		#The z in the self.data_arr is in [cm], so to switch to [m] we use 0.01				
		f = Function()
		if(zSimmetric > 0):
			for iz in range(1,self.zSteps+1):
				[z,r,Ez,Er,E,B] = self.data_arr[iz]
				z = self.Zmin + stepZ*iz
				f.add(-z*0.01,Ez/Ez_max)			
		for iz in range(self.zSteps+1):
			[z,r,Ez,Er,E,B] = self.data_arr[iz]
			z = self.Zmin + stepZ*iz
			f.add(z*0.01,Ez/Ez_max)			
		spline = SplineCH()
		spline.compile(f)
		return spline
