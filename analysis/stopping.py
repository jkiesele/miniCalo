import uproot3 as uproot
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

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

    
def readFile(filename, particles):
    tree = uproot.open(filename)["B4"]
    
    nlayers = np.max(tree["layerNo"].array()[0])+1
    nevents = tree["Nevents"].array()[0]
    energy = tree["true_energy"].array()[0]
    print(nlayers, len(particles))
    
    allarrays=[]
    
    for p in particles:
        branch = tree[p.branchname()+"_stopped"]
        a = branch.array().flatten()
        a = np.asarray(a)
        a = np.reshape(a, [-1,nlayers])
        print(p, np.max(a))
        a = a/float(nevents)
        a = np.sum(a, axis=0) 
        allarrays.append(a)

    return allarrays, nlayers, nevents, energy

def removeCMS(a):
    return a[4:]

particles=[]
with open(G4RHADRONS+"/particles.txt") as f:
    for l in f:
        l = removeEmpty(l.split(' '))
        l[-1] = l[-1][:-1]#remove \n
        p = particle(l[0],l[1],l[-1])
        particles.append(p)
        
        
filename="../../build/out_beta_0.1_0.root"


allarrays,_,_,_=readFile(filename,particles)

for p,a in zip(particles,allarrays):
    plt.plot(a, label=p.name, marker='o',linewidth=None)

plt.yscale("log")
plt.legend()
plt.ylabel("Stopping efficiency")
plt.xlabel("Layer number ( < 4: CMS)")
plt.axvline(x = 3.5, color = 'b') 
plt.show()
plt.close()



for p,a in zip(particles,allarrays):
    s = np.cumsum(removeCMS(a),axis=-1)
    plt.plot(s, label=p.name, marker='o',linewidth=None)
    
    
plt.yscale("log")
plt.ylabel("Cumulative stopping efficiency")
plt.xlabel("Layers passed")
plt.legend()
plt.show()
plt.close()

allars=[]
#out_beta_0.2_0
#now other files

fig, ax1 = plt.subplots()

for b in ["0.025", "0.05", "0.1", "0.2", "0.4", "0.45", "0.47"]:
    arrs,nlayers, nevents, energy = readFile("../../build/out_beta_"+b+"_0.root",particles)
    npa =[]
    for a in arrs:
        print(a.shape,"shape")
        npa.append(np.expand_dims(a,axis=0))
    #this is Npart x layer
    arrs = np.sum(np.concatenate(npa,axis=0),axis=0)#sum over different particle types -> Nlayer
    arrs = np.cumsum(removeCMS(arrs),axis=-1)
    ax1.plot(arrs, label="beta "+b+" E="+str(round(energy,2))+" GeV", marker='o')
    
ax1.set_yscale("log")
ax1.set_ylabel("Cumulative stopping efficiency")
ax1.set_xlabel("Layers passed")
ax1.legend()

ax2 = ax1.twiny()
ax2.set_xlabel("Brass passed [cm]")
ax2.set_xticks( ax1.get_xticks() )
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels([x * 2 for x in ax1.get_xticks()])

plt.show()
fig.savefig("test.png")
#make summary plot



    
    



