import numpy as npy
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os


rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 20
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
floc = "/Users/olson45/Documents/Conference Travel/ISSW32_2019/ISSW32_bjo/figures"


expData = npy.loadtxt('expData2.txt',delimiter=',')

RES = [64,128,256]
syms = {64:'k-.',128:'k--',256:'k-'}

LW = 2.5
for rr in RES:
    simData = npy.loadtxt('simData_%s.txt' % rr)
    half = simData.shape[1]/2
    cp = .5*(npy.abs(simData[0,half:]) + npy.abs(simData[0,half:0:-1]))
    #plt.plot( npy.abs(simData[0,half:]), simData[1,half:], syms[rr],label='$R_c/\Delta x = %s$' % int(rr/10),linewidth=LW)
    #plt.plot( npy.abs(simData[0,:half]), simData[1,:half], syms[rr],label='$R_c/\Delta x = %s$' % int(rr/10),linewidth=LW)
    plt.plot( cp, simData[1,half:], syms[rr],label='$R_c/\Delta x = %s$' % int(rr/10),linewidth=LW)

plt.plot( expData[:,0], expData[:,1], 'ko',label='Exp (Metcalf et al.)',linewidth=LW)

plt.legend()
plt.xlabel(r"$\theta$") 
plt.ylabel(r"$C_p$")
plt.tight_layout(pad=1.0)
plt.savefig(os.path.join(floc,'cyl_exp.pdf'))


plt.show()

