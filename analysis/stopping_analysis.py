import uproot3 as uproot
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import pandas
import glob
import plotly.express as px
import math 


def readFile(filename):
    tree = uproot.open(filename)["B4"]
    
    nlayers = np.max(tree["layerNo"].array()[0])+1
    pbranches = [
    #'mass','beta',
    '1000021_stopped',
    '1000022_stopped',
    '1000993_stopped',
    '1009213_stopped',
    '1009113_stopped',
    '1093122_stopped',
    'm1009213_stopped']
    
    allarrays=[]
    nevents = tree["Nevents"].array()[0]
    beta = tree["beta"].array()[0]
    mass = tree["mass"].array()[0]
    
    for p in pbranches:
        branch = tree[p]
        a = branch.array().flatten()
        a = np.asarray(a)
        a = np.reshape(a, [-1,nlayers])
        a = a/float(nevents)
        a = np.sum(a, axis=0,keepdims=True) 
        allarrays.append(a)

    #print(allarrays)
    #check if any stopped
    a = np.concatenate(allarrays, axis=0)
    #print(a)
    #sum all particle types
    sums = np.sum(a,axis=0)
    sums = sums[4:] #remove CMS
    dsum = sums
    sums = np.cumsum(sums,axis=-1)
    #print(dsum)
    return beta, mass, sums, dsum

def createDF():
    files = glob.glob("out_m*_beta_*.root")
    d={}
    d['mass']=[]
    d['totalabs']=[]
    d['beta']=[]
    d['maxlayer']=[]
    d['energy']=[]
    d['totalabsathalf']=[]
    d['totalabsatquart']=[]
    for f in files:
        try:
            beta,mass,sums,dsum = readFile(f)
        
            #print(dsum)
            maxatlayer = np.argmax(dsum,axis=-1)
            #print(maxatlayer)
            
        
            d['mass'].append(mass/1000.)
            d['totalabs'].append(sums[-1])
            halfidx = sums.shape[-1]//2-1
            #print(halfidx)
            d['totalabsathalf'].append(sums[halfidx])
            d['totalabsatquart'].append(sums[halfidx//2])
            d['beta'].append(beta)
            d['maxlayer'].append(float(maxatlayer))
            gamma = 1 / math.sqrt(1 - beta*beta)
            energy = (gamma-1.)*mass/1000.
            d['energy'].append(energy)
            #d['']
            #plt.plot(range(len(sums)),sums, label='beta='+str(beta))
        except:
            pass
        #break#exit()
        
    
    df = pandas.DataFrame(d,columns=['mass','totalabs','beta','maxlayer','energy','totalabsathalf','totalabsatquart'])    
    df.to_pickle("dataframe.pkl")
    exit()
    

#createDF()


df = pandas.read_pickle("dataframe.pkl")

#exit()

df['logE'] = np.log10(df['energy'])

fig = px.scatter_matrix(df, 
                        dimensions=['mass','totalabs','totalabsathalf','totalabsatquart','beta','maxlayer','energy','logE'], 
                        color="totalabs",
                        #template='plotly_dark'
                        )

fig.write_html("/eos/home-j/jkiesele/www/files/test.html")


#now some standard plots
masses = np.unique(df['mass'])

def makeAbsPlot(fieldname,xfieldname):
    plt.close()
    betas_m, abs_m= {},{}
    for m in masses:
        sel = df.loc[df['mass'] == m]
        betas = np.array(sel[xfieldname])
        abs = np.array(sel[fieldname])
        sorting = np.argsort(betas)
        betas = betas[sorting]
        abs = abs[sorting]
        betas_m[m] = betas
        abs_m[m] = abs
        plt.plot(betas, abs,label=str(int(round(m-0.7)))+' GeV')
    
    plt.legend(title=r'$\tilde{g}$ mass')
    plt.xlabel(r"$\beta$")
    plt.ylabel("Absorption efficiency")
    return betas_m, abs_m
    
   
for x in ['beta','energy']: 
    betas, abs = makeAbsPlot('totalabs',x)
    plt.title("2 m brass")
    plt.ylim([0,0.13])
    if x == 'energy':
        plt.xlabel(r"$E_{kin} = p_{T}$ [GeV]")
        plt.xlim([0,100])
    plt.savefig("totalabs_"+x+".pdf")
    for m in betas.keys():
        np.savetxt("totalabs_"+x+'_'+str(m)+"GeV.txt",betas[m]  ,fmt='%.18e,',header='[',footer=']',comments='')
        np.savetxt("totalabs_"+x+'_'+str(m)+"GeV_abs.txt",abs[m],fmt='%.18e,',header='[',footer=']',comments='')
    continue
    
    
    betas, abs = makeAbsPlot('totalabsathalf',x)
    plt.title("1 m brass")
    plt.ylim([0,0.13])
    if x == 'energy':
        plt.xlabel(r"$E_{kin} = p_{T}$ [GeV]")
        plt.xlim([0,100])
    plt.savefig("totalabsathalf_"+x+".pdf")
    
    betas, abs = makeAbsPlot('totalabsatquart',x)
    plt.title("0.5 m brass")
    plt.ylim([0,0.13])
    if x == 'energy':
        plt.xlabel(r"$E_{kin} = p_{T}$ [GeV]")
        plt.xlim([0,100])
    plt.savefig("totalabsatquart_"+x+".pdf")
    





