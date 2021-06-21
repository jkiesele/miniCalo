

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl



class rod(object):
    def __init__(self, U, pos):
        assert len(pos)>1 and len(pos) < 3
        self.pos = pos
        self.U = U
        


class EField(object):
    def __init__(self,rods, eps_r=1):
        assert len(rods)>0
        self.rods=rods
        self.eps_r = eps_r
    
    def UtoQ(self,myrod):
        return -myrod.U
        
    def getFieldVectors(self, x, y):#x and y np arrays
        
        x = np.expand_dims(x,axis=-1)
        
        y = np.expand_dims(y,axis=-1)

        probepos = np.concatenate([x,y],axis=-1)
        print('probepos a',probepos.shape)
        probepos = np.reshape(probepos, [1, probepos.shape[0]*probepos.shape[1],-1])
        

        eps_0 = 8.854e-12
        pi = 3.1415
        
        Qs = np.array([[self.UtoQ(r)] for r in self.rods]) # R x 1
        pos = np.array([ r.pos for r in self.rods]) #R x 2
        Qs = np.expand_dims(Qs, axis=1) # R x 1 x 1
        pos = np.expand_dims(pos, axis=1) # R x 1 x 2
        
        
        
        distances = np.sqrt(np.sum( (probepos-pos)**2 , axis=-1, keepdims=True)) # R x P x 1
        
        print('distances', distances.shape)
        print('probepos', probepos.shape)
        print('pos', pos.shape)
        print('Qs', Qs.shape)
        
        vecdistances = (probepos-pos)
        #take into account rod sizes
        distances = vecdistances/distances**3
        outside_rod = np.expand_dims(np.any(np.abs(vecdistances) >  0.5, axis=-1),axis=-1)
        distances = np.where( outside_rod, 
                              distances,  0.)
        
        E = 1/(4*pi*eps_0 * self.eps_r) * np.sum( Qs* distances  ,axis=0) #sum over rods: P x 2
        print("E",E.shape)
        
        return E[:,0], E[:,1], np.sqrt( np.sum(E**2, axis=-1) ) #Ex Ey
        
    


rod_dist = 2
#create rods:
x, y = np.meshgrid(np.linspace(rod_dist*(-3), rod_dist*3, 7), 
                   np.linspace(rod_dist*(-3), rod_dist*3, 7))
x = np.expand_dims(x,axis=-1)
y = np.expand_dims(y,axis=-1)
rodpos = np.concatenate([x,y],axis=-1)
rodpos = np.reshape(rodpos, [rodpos.shape[0]*rodpos.shape[1],-1])

print('rodpos',rodpos.shape)


#plot the rods


fig, ax = plt.subplots()
# Add the patch to the Axes
for rxy in rodpos:
    rect = patches.Rectangle(rxy-0.5, 1., 1., linewidth=1, edgecolor='b', facecolor='none')
    ax.add_patch(rect)


#altering charges
rods = []
U=500
for i in range(len(rodpos)):
    r = rod(U,rodpos[i])
    U *= -1
    rods.append(r)

efield = EField(rods)


## create probe grid

x, y = np.meshgrid(np.linspace(rod_dist*(-5), rod_dist*5, 71), 
                   np.linspace(rod_dist*(-5), rod_dist*5, 71))


u,v, Eabs = efield.getFieldVectors(x,y)

print('u,v',u.shape,v.shape)
print('x,y',x.shape,y.shape)
print('Eabs',Eabs.shape)

cm = mpl.cm.cividis
# Plotting Vector Field with QUIVER
plt.quiver(x, y, u, v, color=cm((Eabs-np.min(Eabs))/np.max(Eabs-np.min(Eabs))), pivot='middle')
plt.title('Electric field')
plt.xlabel("y [cm]")
plt.ylabel("z [cm]")
# Setting x, y boundary limits
plt.xlim(-8, 8)
plt.ylim(-8, 8)
  
# Show plot with gird
plt.grid()
plt.show()