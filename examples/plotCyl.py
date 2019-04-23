import numpy as npy
import matplotlib.pyplot as plt


expData = npy.loadtxt('expData.txt',delimiter=',')

RES = [64,128,256]
syms = {64:'k-.',128:'k--',256:'k-'}

LW = 2.5
for rr in RES:
    simData = npy.loadtxt('simData_%s.txt' % rr)
    half = simData.shape[1]/2
    plt.plot( npy.abs(simData[0,half:]), simData[1,half:], syms[rr],label='%s' % rr,linewidth=LW)

plt.plot( expData[:,0], expData[:,1], 'ko',label='Exp',linewidth=LW)


plt.legend()

plt.show()

