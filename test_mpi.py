
from numpy import *
import sys
import os
import time
from init_functions import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
nprocs = comm.size
myrank = comm.rank


x = asarray(myrank)
total = asarray(0)
comm.Allreduce(x, total, op=MPI.SUM)

print 'total is ',total
