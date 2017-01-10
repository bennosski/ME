
from mpi4py import MPI
from numpy import *

comm = MPI.COMM_WORLD
nprocs = comm.size
myrank = comm.rank

print "Hello! I'm rank %d from %d running in total..." % (comm.rank, comm.size)

mysum = array([0.])

#loop from myrank by nprocs until loop is done
for i in range(myrank,13,nprocs):
    print "printing i ",i," from rank ",comm.rank
    mysum += array([i])
    
total = array([0.])
comm.Allreduce(mysum, total, op=MPI.SUM)

if myrank==0:
    print "total is ",total
                    
