import subprocess
import shutil
import threading
import signal
import datetime
import os
import copy
import math
import sys
from operator import attrgetter

print os.path.abspath(os.path.join(os.getcwd(), '..', '..', '..', 'dataparser'))
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..', '..', '..', 'dataparser')))

print dir()

from dataparser import parse_file

class InputData(object):
    global_solver = 'BJ'
    local_solver = 'vbarms'
    nlev = 3
    droptol = 1e-3
    lfil = 100000
    tau = 0.8
    maxits = 1000
    cosine = 1

    def __str__(self):
        return """%s           parallel preconditioner (RAS, SCHUR, BJ)
%s          local preconditioner (ilu0, iluk, ilut, arms, vbarms, vbilut, vbiluk)
0.01            eps(tolerance for inner iteration)
1.0e-6          eps1(tolerance for outer iteration)
%d               nlev(number of levels for ARMS, ILUK and VBARMS, VBILUK)
0               non-symmetric permutation - ddPQ or not (1 = ddPQ)
250             bsize(block size for block independent sets)
0.40            tolind(tolerance used in independent set)
500              im (krylov subspace size for outer iteration)
%d            maxits (outer fgmres)
10              kim (krylov subspace size for inner iteration)
10              itsgmr(inner fgmres)
0               nonsymmetric permutations for interlevel blocks (1:yes, 0:no)
0               permutations of columns for interlevel blocks (1:yes, 0:no)
1               row scaling for interlevel blocks (1:yes, 0:no)
1               column scaling for interlevel blocks (1:yes, 0:no)
0               nonsymmetric perm.s for last Schur complement (1:yes, 0:no)
0               permutations of columns for last Schur complement (1:yes, 0:no)
1               row scaling for last Schur complement (1:yes, 0:no)
1               column scaling for last Schur complement (1:yes, 0:no)
%d    lfil0(ilut, iluk and arms for lfil[0-3])
%d    lfil4(schur)
%d    lfil5(ILUT L, ILUT U)	
%.3e          droptol0(droptol[0=3], L, U, L^{-1}F, EU^{-1}
%.3e          droptol4(for schur complements at each level)
%.3e          droptol5(for ILUT in last level Schur Complement)
%f            eps_init(For block version local solver(vbarms and vbilut) only, angle tolerance of init-block subroutine)
%d            cosine
##

""" % (self.global_solver, self.local_solver, self.nlev, self.maxits, self.lfil, self.lfil, self.lfil, self.droptol, self.droptol, self.droptol, self.tau, self.cosine)

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            print 'Thread started'
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            print 'Thread finished at ', str(datetime.datetime.now())

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print 'Terminating process'
            #~ self.process.send_signal(signal.CTRL_C_EVENT)
            subprocess.call('killall '+self.cmd.partition(' ')[0], shell=True)
            thread.join()
        #~ timer = threading.Condition()
        #~ timer.acquire()
        #~ timer.wait(0.1)
        print self.process.returncode
        return self.process.returncode
#~ 
#~ class Alarm(Exception):
    #~ pass
#~ 
#~ def alarm_handler(signum, frame):
    #~ raise Alarm
#~ 
#~ class Command(object):
    #~ def __init__(self, cmd):
        #~ self.cmd = cmd
#~ 
    #~ def run(self, timeout):
        #~ proc = subprocess.Popen(self.cmd, shell=True)
        #~ signal.signal(signal.SIGALRM, alarm_handler)
        #~ signal.alarm(timeout)
        #~ try:
            #~ stdoutdata, stderrdata = proc.communicate()
            #~ signal.alarm(0)
        #~ except Alarm:
            #~ print "Oops, taking too long!"
            #~ proc.kill()

def write_input(data):
    f = open('inputs', 'w')
    f.write(str(data))
    #~ f.write('bj              parallel preconditioner (RAS, SCHUR, BJ)\n'+method+\
    #~ '          local preconditioner (ilu0, iluk, ilut, arms, vbarms, vbilut, vbiluk)\n'+inputs)
    f.close()

def get_fname(matrix, data, nproc):
    return 'output/' + matrix.split('/')[4].split('.')[0] + '_' + data.local_solver +  '_' + str(nproc) + '.txt'

def write_debug(text):
    debugfile = open('debug.txt', 'a')
    print text
    debugfile.write(text)
    debugfile.close()

def test_method(matrix, data, nproc=1):
    write_input(data)

    if not data.local_solver.startswith('vb'):
        exe = './dd-mtx-dse.ex'
    else:
        exe = './Bdd-mtx-dse.ex'

    #~ if nprocs = None:
        #~ nprocs = [1]
    #~ for n in nprocs:
    n = nproc

    fname = get_fname(matrix, data, n)
    c = Command('mpiexec -n ' + str(n) + ' ' + exe)
    ret = c.run(600)

    #~ shutil.copy('output.txt', fname)
    f = open('output.txt', 'r')
    l = f.readlines()
    f.close()
    f = open(fname, 'w')
    f.write(''.join(l[-13:]))
    f.close()
    new_data = parse_file(fname)

    if ret != 0:
        new_data.isgood = False

    pp = matrix.split('/')[4].split('.')[0] + ' with '+ method + ' and ' + str(data.droptol)+ ' and '+ str(int(data.nlev))+ ', t='+ str(new_data.total_time) + ', mem=' + str(new_data.memory) + ', tau='+ str(data.tau)+', ret=' + str(ret) +'\n'
    write_debug(pp)

    return new_data

def get_mem(fc):
    return fc.memory if not isinstance(fc.memory, basestring) else 1e-100

if __name__ == '__main__':
    mfile = open('matrices.txt', 'r')
    matrices = mfile.readlines();
    mfile.close()

    debugfile = open('debug.txt', 'w')
    debugfile.write('')
    debugfile.close()

    nproc = 1
    overwrite = True

    #~ f = open('inputs_short', 'r')
    #~ inputs = f.read()
    #~ f.close()

    #f = open('results.txt', 'w')
    #f.write('')
    #f.close()
    for m in matrices:
        f = open('matfileReal', 'w')
        f.write(m+'\n##')
        f.close()

        fdict = {}
        fvdict = {}
        data = InputData()

        for method in ['vbarms']:
            data.local_solver = method
            fname = get_fname(m, data, nproc)
            if os.path.isfile(fname) and not overwrite:
                fv = parse_file(fname)
                continue

            for tau in (0.6, 0.7, 0.8, 0.9, 1.0):
                data.tau = tau

                data.droptol = 1
                fa = test_method(m, data)
                fdict[fa] = copy.copy(data)

                data.droptol = 1e-8
                fb = test_method(m, data)
                fdict[fb] = copy.copy(data)

                testing = fa.isgood or fb.isgood
                while testing:
                    data.droptol = math.pow(10, (math.log10(fdict[fb].droptol) + math.log10(fdict[fa].droptol)) / 2)
                    fc = test_method(m, data)
                    fdict[fc] = copy.copy(data)

                    if not fa.isgood or (fc.isgood and fc.total_time < fa.total_time):
                        fa = fc
                    else:
                        fb = fc
                    testing = (fa.isgood or fb.isgood) and (math.log10(fdict[fb].droptol) - math.log10(fdict[fa].droptol) <= -1)
                fl = sorted([fa, fb ,fc], key=attrgetter('total_time'))
                fvdict[tau] = fl[0]

            fv = sorted(fvdict.values(), key=attrgetter('total_time'))[0]
            f = open(os.path.join('input', 'inputs_'+m.split('/')[4].split('.')[0]+'_'+method), 'w')
            f.write(str(fdict[fv]))
            f.close()

            data = copy.copy(fdict[fv])
            data.maxits = 1000
            test_method(m, data)

        for method in ['vbarmsold', 'ilut', 'iluk', 'arms', 'vbilut', 'vbiluk']:
            data.local_solver = method
            fname = get_fname(m, data, nproc)
            if os.path.isfile(fname) and not overwrite:
                continue

            data.nlev = 3
            data.maxits = 1

            if method.endswith('k'):
                data.nlev = 0
            else:
                data.droptol = 1
            fa = test_method(m, data)
            fdict[fa] = copy.copy(data)
            if method.endswith('k'):
                data.nlev = 8
            else:
                data.droptol = 1e-8
            fb = test_method(m, data)
            fdict[fb] = copy.copy(data)
            testing = fa.isgood or fb.isgood
            testing = testing and (not fa.isgood or fa.memory < fv.memory)
            testing = testing and (not fb.isgood or fb.memory > fv.memory)

            while testing:
                if method.endswith('k'):
                    data.nlev = math.floor((fdict[fb].nlev + fdict[fa].nlev) / 2)
                else:
                    data.droptol = math.pow(10, (math.log10(fdict[fb].droptol) + math.log10(fdict[fa].droptol)) / 2)
                fc = test_method(m, data)
                fdict[fc] = copy.copy(data)

                if not fa.isgood and not fc.isgood:
                    fa = fc
                elif fc.isgood and get_mem(fc) < get_mem(fv):
                    fa = fc
                else:
                    fb = fc
                testing = (fa.isgood or fb.isgood) and abs(math.log10(get_mem(fv)) - math.log10(get_mem(fc))) > 0.1
                if method.endswith('k'):
                    testing = testing and (fdict[fb].nlev - fdict[fa].nlev > 1)
                #~ else:
                    #~ testing = (fa.isgood or fb.isgood) and (math.log10(fb_data.droptol) - math.log10(fa_data.droptol) <= -0.1)

            if abs(math.log10(get_mem(fv)) - math.log10(get_mem(fa))) < abs(math.log10(get_mem(fv)) - math.log10(get_mem(fb))):
                fc = fa
            else:
                fc = fb

            f = open(os.path.join('input', 'inputs_'+m.split('/')[4].split('.')[0]+'_'+method), 'w')
            f.write(str(fdict[fc]))
            f.close()

            data = copy.copy(fdict[fc])
            data.maxits = 1000
            test_method(m, data)
