
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



L = np.array([3, 4, 6, 8, 10, 12],dtype='float')
L -= 2. #remove rods
L /= 200 #individual spacing
L *= 100. #in cm

relLAr = np.array([0.22, 0.35, 0.52, 0.69, 0.80, 0.85])
relLArDelta =     [0.14, 0.1, 0.05, 0.05, 0.04, 0.03]

#relLArTrue = np.array([0.29,0.5,0.62,0.72,0.77])
#relLArTrueDelta = [0.08,0.13,0.12,0.1,0.13]

plt.errorbar(L, relLAr, yerr=relLArDelta,marker='o',linewidth=0,elinewidth=2,lineStyle=None)
plt.xlabel("rod spacing [cm]")
plt.ylabel(r"$E_{LAr}/(E_{LAr}+E_{rod})$")
plt.savefig("relEdep.pdf")
plt.close()

exit()

plt.errorbar(L, relLArTrue, yerr=relLArTrueDelta,marker='o',linewidth=0,elinewidth=2,lineStyle=None)
plt.xlabel("rod spacing [cm]")
plt.ylabel(r"$E_{LAr}/E_{hadron}$")
plt.savefig("relEinLAr.pdf")
plt.close()


