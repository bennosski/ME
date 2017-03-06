import sys, os
from numpy import *

#g = float(sys.argv[1])
#omega = 0.1
#N = 64
#mu = -N*g**2/omega**2

import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

mus = arange(-4.0,0.0,0.25)
mus = [-4.0,-3.75,-3.5]

ct = 0
cmd = "echo hi"
for mu in mus:    

     print "setting up input file with mu ",mu

     # sed "/chemical/c\mu = $2  #chemical potential" sample_input_forward.txt > "# input_mu_search_$1.txt'
     cmd += "; cp input temp1"
     cmd += """; sed "/mu/c\mu     # %1.15f"""%mu + """ " temp1 > input%d"""%ct
     cmd += "; sbatch submit input%d"%ct

     ct += 1

bash_command(cmd)
