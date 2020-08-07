
import numpy as np
import uproot
import ROOT
import glob

path="/eos/home-d/dorigo/Muon_E_loss/sample_3/"

def readFile(filename):

    tree = uproot.open(filename)["B4"]
    true_energy = np.array( list(tree["true_energy"].array()) ,dtype='float32')
    return true_energy



def checkAllFiles():
    files = glob.glob(path+"*.root")
    allgood=[]
    for f in files:
        a = readFile(f)
        un,c = np.unique(a, return_counts=True)
        if len(un[c>1]):
            print (f, un[c>1])
            tree = uproot.open(f)["B4"]
            total_dep_en = np.array( list(tree["total_dep_energy"].array()) ,dtype='float32')
            for dd in un[c>1]:
                print(total_dep_en[a==dd])
        else:
            allgood.append(f)
    return allgood, files
    
allgood,files = checkAllFiles()
print('all files checked', len(allgood)/len(files) * 100,'% ok')
import random
random.shuffle(allgood)

groups = zip(*[iter(allgood)]*2)

for fa, fb in groups:
    a = readFile(fa)
    b = readFile(fb)
    
    all = np.concatenate([a,b],axis=0)
    
    un,c = np.unique(all, return_counts=True)
    if len(un[c>1]):
        print(len(un), len(all),un[c>1])