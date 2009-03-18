import orbit_mpi
from orbit_mpi import *

class printf:

    def __init__(self,f_name,*title):

        self.f_name = f_name
        rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
        if (rank == 0):
            f_out = open(self.f_name,"w")
            for i in range(len(title)):
                f_out.write("%19s"%title[i])
            f_out.write("\n")
            f_out.close()



    def fdata(self,*data):

        rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
        if (rank == 0):
            f_out = open(self.f_name,"a")
            for i in range(len(data)):
                f_out.write("%19.4f"%data[i])
            f_out.write("\n")
            f_out.close()

