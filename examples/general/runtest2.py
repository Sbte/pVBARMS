import subprocess
import shutil

mfile = open('matrices.txt', 'r')
matrices = mfile.readlines();
mfile.close()

f = open('results.txt', 'w')
f.write('')
f.close()
for m in matrices:
    f = open('matfileReal', 'w')
    f.write(m+'\n##')
    f.close()

    shutil.copyfile('inputs_vbarms', 'inputs')
    subprocess.call('mpiexec -n 3 ./Bdd-mtx-dse.ex >> results.txt', shell=True)
#~ 
    #~ shutil.copyfile('inputs_arms', 'inputs')
    #~ subprocess.call('mpiexec -n 3 ./dd-mtx-dse.ex >> results.txt', shell=True)
