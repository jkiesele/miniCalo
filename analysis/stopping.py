import uproot3 as uproot
import numpy as np
import matplotlib.pyplot as plt
import os

G4RHADRONS = os.getenv("G4RHADRONS")

#read possible particles
def removeEmpty(al):
    alo=[]
    for l in al:
        if len(l):
            alo.append(l)
    return alo

class particle(object):
    
    def __init__(self,pdgid=0,mass=0,name=""):
        self.mass =  float(mass)
        self.pdgid = int(pdgid)
        self.name = name
        
    def branchname(self):
        if self.pdgid<0:
            s = str(self.pdgid)
            s="m"+s[1:]
            return s
        else:
            return str(self.pdgid)
    
    def __str__(self, *args, **kwargs):
        return self.name + ", " +str(self.pdgid) +", mass: "+ str(self.mass)

    

particles=[]
with open(G4RHADRONS+"/particles.txt") as f:
    for l in f:
        l = removeEmpty(l.split(' '))
        l[-1] = l[-1][:-1]#remove \n
        p = particle(l[0],l[1],l[-1])
        particles.append(p)
        
        
filename="../../build/out0.root"

tree = uproot.open(filename)["B4"]

nlayers = np.max(tree["layerNo"].array()[0])+1
print(nlayers, len(particles))

for p in particles:
    branch = tree[p.branchname()+"_stopped"]
    a = branch.array().flatten()
    a = np.asarray(a)
    a = np.reshape(a, [-1,nlayers])
    print(a.shape)
    a = a/float(a.shape[0])
    a = np.sum(a, axis=0) 
    print (a.shape)
    plt.plot(a, label=p.name, marker='o',linewidth=None)

plt.yscale("log")
plt.legend()
plt.show()
    
    



