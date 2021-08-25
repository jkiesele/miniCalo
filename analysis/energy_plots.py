
import numpy as np
import matplotlib.pyplot as plt


import matplotlib
import matplotlib.pyplot as plt
import math
fontsize=12

axes = {'labelsize': fontsize,
        'titlesize': fontsize}

matplotlib.rc('axes', **axes)
matplotlib.rc('legend',fontsize=fontsize)
matplotlib.rc('xtick',labelsize=fontsize)
matplotlib.rc('ytick',labelsize=fontsize)



L = np.array([4, 6, 8, 10, 12],dtype='float')
L -= 2. #remove rods
L /= 200 #individual spacing
L *= 100. #in cm

relLAr = np.array([0.324, 0.56, 0.68, 0.77, 0.83])
relLArDelta = [0.06, 0.1, 0.06, 0.05, 0.04]

relLArTrue = np.array([0.29,0.5,0.62,0.72,0.77])
relLArTrueDelat = [0.08,0.13,0.12,0.1,0.13]

plt.errorbar(L, relLAr, yerr=relLArDelta,marker='o',linewidth=0,elinewidth=2,lineStyle=None)
plt.xlabel("rod spacing [cm]")
plt.ylabel(r"$E_{LAr}/(E_{LAr}+E_{rod})$")
plt.savefig("relEdep.pdf")
plt.close()

plt.errorbar(L, relLArTrue, yerr=relLArTrueDelat,marker='o',linewidth=0,elinewidth=2,lineStyle=None)
plt.xlabel("rod spacing [cm]")
plt.ylabel(r"$E_{LAr}/E_{hadron}$")
plt.savefig("relEinLAr.pdf")
plt.close()


