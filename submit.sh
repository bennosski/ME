#!/bin/bash
#SBATCH --job-name=ME
#SBATCH --time=00:10:00
#SBATCH --constraint="2.60GHz"
#SBATCH --partition=iric
#SBATCH --qos=iric
#SBATCH --ntasks=192
#SBATCH --cpus-per-task=1

##after doing a pip install --user mpi4py
module load python/2.7.5
export PYTHONPATH=~/.local/lib/python2.6/site-packages/:$PYTHONPATH

echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
echo "^^^^^^^start: `date`"
srun --ntasks=192 --cpus-per-task=1 python ME1_mpi.py>log
echo "^^^^^^^end: `date`"
