#!/bin/bash
#SBATCH --job-name=ME
#SBATCH --time=00:05:00
#SBATCH --constraint="2.60GHz"
#SBATCH --partition=iric
#SBATCH --qos=iric
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1

##after doing a pip install --user mpi4py
module load python/2.7.5
export PYTHONPATH=~/.local/lib/python2.6/site-packages/:$PYTHONPATH

echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
echo "^^^^^^^start: `date`"
srun --ntasks=20 --cpus-per-task=1 python ME_forward_mpi.py>log
echo "^^^^^^^end: `date`"
