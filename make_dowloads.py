import sys
import os
import subprocess
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])

dirpath = sys.argv[1]
folders = os.listdir(dirpath)

for folder in folders:
    if not os.path.isdir(dirpath+folder):
        continue

    if "_d" in folder:
        continue

    print folder
    
    folder_d = folder + "_d/"
    
    if not os.path.exists(dirpath+folder_d):
        os.makedirs(dirpath+folder_d)

    cmd  = "cp "+dirpath+folder+"/dens "+dirpath+folder_d
    cmd += "; cp "+dirpath+folder+"/xsc* "+dirpath+folder_d
        
    bash_command(cmd)
    
    

