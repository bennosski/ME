
########
HOLSTEIN
########

ME_holstein_mpi.py    :::: solves imaginary frequency holstein problem
Marsiglio_holstein.py :::: uses output of ME_holstein_mpi.py in /data_holstein folder and does Marsiglio iterations


#######

THESE ARE CHECKED

ME_forward.py            ::: solves imaginary frequency forward scattering problem
ME_forward_mpi.py       ::: uses MPI. Is the latestest used verions
Marsiglio_forward.py  ::: uses output of ME1_mpi.py in /data folder and does Marsiglio iterations



setup notes:
      ml python/2.7.5
      how to load MPI on sherlock:
      export PYTHONPATH=~/.local/lib/python2.6/site-packages/:$PYTHONPATH
