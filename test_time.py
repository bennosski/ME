import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
nprocs = comm.size
myrank = comm.rank


s1 = time.time()


for i in range(myrank, 20, nprocs):
    #print "hi"
    time.sleep(1.0)
    #print "end"

print myrank, ' took ',time.time()-s1
