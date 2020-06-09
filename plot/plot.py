from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LightSource

import matplotlib as mpl
import matplotlib.cm as cm
import uproot



def cuboid_data(pos, size=(1,1,3)):
    # code taken from
    # https://stackoverflow.com/a/35978146/4124317
    # suppose axis direction: x: to left; y: to inside; z: to upper
    # get the (left, outside, bottom) point
    o = [a - b / 2 for a, b in zip(pos, size)]
    # get the length, width, and height
    l, w, h = size
    x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],  
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],  
         [o[0], o[0] + l, o[0] + l, o[0], o[0]],  
         [o[0], o[0] + l, o[0] + l, o[0], o[0]]]  
    y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],  
         [o[1], o[1], o[1] + w, o[1] + w, o[1]],  
         [o[1], o[1], o[1], o[1], o[1]],          
         [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]   
    z = [[o[2], o[2], o[2], o[2], o[2]],                       
         [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],   
         [o[2], o[2], o[2] + h, o[2] + h, o[2]],               
         [o[2], o[2], o[2] + h, o[2] + h, o[2]]]               
    return np.array(x), np.array(y), np.array(z)

def plotCubeAt(pos=(0,0,0),size=(1,1,3),color=0,ax=None, alpha=1, m=None):
    # Plotting a cube element at position pos
    if ax !=None:
        X, Y, Z = cuboid_data( pos, size )
        ax.plot_surface(X, Z, Y, color=m.to_rgba(color), rstride=1, cstride=1, alpha=alpha,
                        antialiased=True, shade=False)


#ax.set_aspect('equal')
getdata=True

def readandshape(tree,varname):
    entry = np.array( list(tree[varname].array()) ,dtype='float32')
    return np.reshape(entry, [-1, 30,30,60,1])


tree = uproot.open("out.root")["B4"]

rechit_energy = readandshape(tree, "rechit_energy")
rechit_x   = readandshape(tree, "rechit_x")
rechit_y   = readandshape(tree, "rechit_y")
rechit_z   = readandshape(tree, "rechit_z")
rechit_vxy =readandshape(tree, "rechit_vxy")
rechit_vz  = readandshape(tree, "rechit_vz")

true_energy = np.array( list(tree["true_energy"].array()) ,dtype='float32')
#print(true_energy.shape)
#exit()

calo = np.concatenate([rechit_energy,
                       rechit_x   ,
                       rechit_y   ,
                       rechit_z   ,
                       rechit_vxy ,
                       rechit_vz  ],axis=-1)
# e, x, y, z
def plotevent(event, arr, ax, iscalo=True):
    usearr = arr[event]
    
    scaled_emax = np.log(np.max(usearr[:,:,0])+1)
    norm = mpl.colors.Normalize(vmin=0, vmax=scaled_emax)
    cmap = cm.hot
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    alpha = 0.5
    if not iscalo:
        alpha=0.01
    for i in range(usearr.shape[0]):
        for j in range(usearr.shape[1]):
            for k in range(usearr.shape[2]):
                e = np.log(usearr[i,j,k,0]+1)
                x = usearr[i,j,k,1]
                y = usearr[i,j,k,2]
                z = usearr[i,j,k,3]
                dxy_indiv = usearr[i,j,k,4]
                dz = usearr[i,j,k,5]
                alpha = (e+0.005)/(scaled_emax+0.005)
                if alpha<0.0005:
                    continue
                plotCubeAt(pos=(x,y,dz/2.+z), size=(dxy_indiv,dxy_indiv,dz), color=e, ax=ax, m=m, alpha=alpha)
    
    
for event in range(10):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel("x [mm]")
    ax.set_zlabel("y [mm]")
    ax.set_ylabel("z [mm]")
    ax.grid(False)
    
    en = true_energy[event]
    ax.text2D(0.05, 0.95, "Energy="+str(en)+" GeV", transform=ax.transAxes)
    
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    print('plotting...')
    plotevent(event,calo,ax,True)
    plt.tight_layout()
    print('show')
    #plt.show()
    plt.savefig("cpion_"+str(event)+".png")


