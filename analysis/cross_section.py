import math

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from muon_spectrum import getMuonsPerSecAboveThresholdBuffered
import styles
import os
os.system('mkdir -p plots')

def hasDecayedProb(t, lifetime):
    scaledlt = t*math.log(2)
    return 1. - np.exp( -t/lifetime )

class massClass(object):
    def __init__(self, mass, position : int=0):

        m = [5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560]
        theoryCrossSectionInFb = [
            2.897e12,
            3.473e11,
            3.236e10,
            2.422e9,
            1.443e8,
            6.615e6,
            2.163e5,
            4.221e3,
            3.254e1,
            2.863e-2,
        ]

        #acc times eff taken from plotHist.cc in Jasmine's code
        if position==0:
            efficiency = [
                7.21837e-05,
                0.0002612,
                0.000513178,
                0.000759026,
                0.000874342,
                0.00104653,
                0.000826617,
                0.000781483,
                0.000653309,
                0.000540727,
            ]
            #linear interpolation between points
            self.totalefficiency = np.interp(mass, m, efficiency)

        elif position==1:
            efficiency = [
                0.00573284,
                0.00644855,
                0.00423486,
                0.00293476,
                0.00196878,
                0.00161342,
                0.000918323,
                0.000779115,
                0.00046486,
                0.000399809,
            ]
            #linear interpolation between points
            self.totalefficiency = np.exp(np.interp(mass, m, np.log(efficiency)))


        self.crosssection = np.exp(np.interp(mass, m, np.log(theoryCrossSectionInFb)))
        self.mass = mass            

    def eff(self):
        return self.totalefficiency

    def xs(self):
        return self.crosssection

    def m(self):
        return self.mass
    
def lumiInTimeFb(
        t_seconds,
        effective_collision_time = 0.5, #half of the time collisions
        L_int :float = 2.5e34 #on average half of peak lumi
        ):
    
    return t_seconds*effective_collision_time*L_int*1e-39
    

#totalHLXCheck = lumiInTimeFb(3600*24*365*9)
#as a check this gives approx 3-4ab-1, given 13 years of operation and assuming 4 years no collisions due to stops
#print('one day of data', lumiInTimeFb(3600*24)) #approx 1fb-1

#HL-LHC should be O(1fb)/day very roughly
def crossSectionNeededFb( 
        mass : massClass,
        lifetime_seconds : float,
        target_events_to_be_detected :float,
        absorption_time_days,
        construction_time_days,
        detection_time_days,
        ):
    
    abstime_seconds = absorption_time_days*3600.*24.
    
    absdecay_seconds = absorption_time_days*3600.*24. / 2.#average
    consdecays_seconds = construction_time_days*3600.*24.
    detdecays_seconds = detection_time_days*3600.*24.
    
    decays_before_detection_prob = hasDecayedProb(absdecay_seconds+consdecays_seconds,lifetime_seconds)
    decays_total_prob = hasDecayedProb(absdecay_seconds+consdecays_seconds+detdecays_seconds,lifetime_seconds)
    
    decays_during_detection_prob = decays_total_prob-decays_before_detection_prob
    
    events_needed = target_events_to_be_detected / (decays_during_detection_prob * mass.eff())
    
    return events_needed/lumiInTimeFb(abstime_seconds)
    

def nEventsDetectable( 
        mass : massClass,
        lifetime_seconds : float,
        absorption_time_days,
        construction_time_days,
        detection_time_days,
        ):
    
    abstime_seconds = absorption_time_days*3600.*24.
    
    absdecay_seconds = absorption_time_days*3600.*24. / 2.#average
    consdecays_seconds = construction_time_days*3600.*24.
    detdecays_seconds = detection_time_days*3600.*24.
    
    decays_before_detection_prob = hasDecayedProb(absdecay_seconds+consdecays_seconds,lifetime_seconds)
    decays_total_prob = hasDecayedProb(absdecay_seconds+consdecays_seconds+detdecays_seconds,lifetime_seconds)
    
    decays_during_detection_prob = decays_total_prob-decays_before_detection_prob
    #print("for lifetime "+str(lifetime_seconds/3600/24)+" days, decays_before_detection_prob is: "+str(decays_before_detection_prob)+", decays_during_detection_prob is: "+str(decays_during_detection_prob)+", decays_total_prob is: "+str(decays_total_prob))

    #N = acc*xs*Lint
    N = decays_during_detection_prob * mass.eff() * mass.xs() * lumiInTimeFb(abstime_seconds)
    #print("For mass "+str(mass.m())+" GeV, the number of events is: "+str(N))
    return N



########### bkg calculations
def comb(n,k):
    import math
    return math.factorial(n) // math.factorial(k) // math.factorial(n - k)


#cache at least M out of N
at_least_k_of_n_cache={}
def at_least_k_of_n(k, n, eff):
    import math
    global at_least_k_of_n_cache
    if (k,n) in at_least_k_of_n_cache.keys():
        return at_least_k_of_n_cache[(k,n)]
    
    def k_of_n(kk,nn,efff):
        return comb(nn, kk)*efff**kk*(1-efff)**(nn-kk)
    out=0
    for i in range(k,n+1):
        out += k_of_n(i,n,eff)
    at_least_k_of_n_cache[(k,n)]=out
    return out
    

# number of cosmic background events, given certain detection time
def calcNBkgEvents(
        detection_time_days,
        cosmicRate = 1, # events/sec/m^2
        effRpc = 0.99, # efficiency of 1 RPC
        nRpcs = 4, #number of RPCs
        M = 2, #number of RPCs that detected the muon, must be at least 2
        absorberSurfaceArea = 4, #m^2 (2m x 2m)
        sameOnTopAndBottom=True
        ):

    detection_time_seconds = detection_time_days*3600.*24.
    N = nRpcs - M #number of RPCs that didn't detect the muon
    
    mod = 1. - at_least_k_of_n(M,nRpcs,effRpc)
    if sameOnTopAndBottom:
        if M%2 or N%2:
            raise ValueError("can't calculate odd number of RPCs with sameOnTopAndBottom=True")
        #at least M//2 out of nRpcs//2 on the top AND on the bottom
        mod = 1. - (at_least_k_of_n(M//2,nRpcs//2,effRpc) * at_least_k_of_n(M//2,nRpcs//2,effRpc))


    return cosmicRate * mod * absorberSurfaceArea * detection_time_seconds

########### THE STUFF BELOW IS JUST FOR TESTING


def makeplots(position):
    
    #efficiency vs mass
    m = [5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560]
    if position==0:
        e = [
            7.21837e-05,
            0.0002612,
            0.000513178,
            0.000759026,
            0.000874342,
            0.00104653,
            0.000826617,
            0.000781483,
            0.000653309,
            0.000540727,
        ]
        massOfAveEff = 1280
    elif position==1:
        e = [
            0.00573284,
            0.00644855,
            0.00423486,
            0.00293476,
            0.00196878,
            0.00161342,
            0.000918323,
            0.000779115,
            0.00046486,
            0.000399809,
        ]
        massOfAveEff = 40
        
    figEff, axEff = plt.subplots()
    eff = axEff.plot(m, e, 'bo')
            
    axEff.set_xlabel("Mass [GeV]")
    axEff.set_ylabel("Acceptance times efficiency")
                        
    plt.savefig('plots/efficiencyVsMass_pos'+str(position)+'.pdf')
    
    #for lifetime_days in [7, 30, 365]:
    for lifetime_days in [7]:
        #for NRpcsFiring in [2, 4]:
        for NRpcsFiring in [2]:
            for construction_time_days in [7]:
        
                x, y = np.meshgrid(np.linspace(lifetime_days/7., 15*lifetime_days, 200), 
                                   np.linspace(lifetime_days/7., 30.*lifetime_days, 200))
            
                z = crossSectionNeededFb(
                    massClass(massOfAveEff,position),
                    lifetime_seconds = 3600.*24*lifetime_days,
                    target_events_to_be_detected = 3.,
                    absorption_time_days = x,
                    construction_time_days = construction_time_days,
                    detection_time_days = y
                )
            
                fig, ax = plt.subplots()
                minzidx = np.unravel_index(np.argmin(z),shape=z.shape) 
            
                print('min',x.shape,)
            
                z = np.log10(z)            
            
                c = ax.pcolormesh(x, y, z, cmap='viridis', vmin=np.min(z), vmax=np.max(z))
            
                # set the limits of the plot to the limits of the data
                ax.axis([x.min(), x.max(), y.min(), y.max()])
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlabel("Absorption time [days]")
                ax.set_ylabel("Detection time [days]")
                ax.set_title('lifetime '+str(lifetime_days)+', construction '+str(construction_time_days)+' [days]')
                fig.colorbar(c, ax=ax,label='Log10 cross section [fb]')
                ax.scatter(x[minzidx],y[minzidx])
                
            
                plt.savefig("plots/times_"+str(lifetime_days)+'_'+str(construction_time_days)+"_pos"+str(position)+'.pdf')


                #sensitivity plot:
                xPoints, yPoints = np.meshgrid(np.logspace(0,np.log10(3000)),np.logspace(0,np.log10(3000)))
                #print (xPoints, yPoints)
                nEvents = np.zeros((len(xPoints),len(xPoints)))
                nBkgEvents = np.zeros_like(nEvents)

                for j in range(len(xPoints)):
                    for i in range(j+1):
                        if yPoints[i,j]< xPoints[i,j]-3: #assume deltaM is 3 GeV, which is the minimum energy we think we can detect
                            nEvents[i,j] = nEventsDetectable(
                                massClass(xPoints[i,j],position),
                                lifetime_seconds = 3600.*24*lifetime_days,
                                absorption_time_days = 2*lifetime_days,
                                construction_time_days = construction_time_days,
                                detection_time_days = 30
                            )
                            
                            deltam = xPoints[i,j] - yPoints[i,j]
                            nBkgEvents[i,j] = calcNBkgEvents(
                                detection_time_days = 30,
                                M = NRpcsFiring,
                                cosmicRate = getMuonsPerSecAboveThresholdBuffered(
                                    deltam/2.,
                                    1. #for 1 m2 since it's already in the function
                                    )
                                )

                            #print("for gluino mass: "+str(xPoints[i,j])+", nEvents["+str(i)+","+str(j)+"] is: "+str(nEvents[i,j]))
                            #if i>10:
                            #    exit()

                nSigEvents = np.where(nEvents==0, np.nan, nEvents)

                #significance = S/sqrt(S+B)
                signif = nSigEvents/np.sqrt(nSigEvents+nBkgEvents)
            
                fig2, ax2 = plt.subplots()
                        
                c2 = ax2.pcolormesh(xPoints, yPoints, signif, cmap='viridis', vmin=0, vmax=5)
            
                # set the limits of the plot to the limits of the data
                ax2.axis([np.min(xPoints), np.max(xPoints), np.min(yPoints), np.max(yPoints)])
                ax2.set_xlabel(r"$m_{\widetilde{g}}$"+" [GeV]")
                ax2.set_ylabel(r"$m_{\widetilde{\chi}^0}$"+" [GeV]")
                ax2.set_xscale('log')
                ax2.set_yscale('log')
                fig2.colorbar(c2, ax=ax2,label='$S/\sqrt{S+B}$')

                #triangle showing kinematically forbidden area
                tri = Polygon(np.array([[np.min(xPoints),np.min(yPoints)], [np.max(xPoints),np.max(yPoints)], [np.min(xPoints),np.max(yPoints)]]), closed=False, fc='lightgrey', ec=None)
                ax2.add_patch(tri)
                #ax2.text(20, 27, "Kinematically forbidden", rotation=44, rotation_mode='anchor')
                ax2.text(20, 25, "Kinematically forbidden", rotation=43, rotation_mode='anchor')

                #superimpose CMS observed limits (doi:10.1007/JHEP05(2018)127, https://www.hepdata.net/record/ins1645630?version=1&table=Table%2018)
                cmsFile = np.genfromtxt('HEPData-ins1645630-v1-Table_18.csv', delimiter=",", names=["x", "y"])
                cms = plt.plot(cmsFile['x'],cmsFile['y'], color='red')
                ax2.text(370, 1.2, "CMS 95% CL observed limits", color='red', rotation=90, rotation_mode='anchor')

            
                ax2.text(1.2, 200,
                         r'$\widetilde{g} \to g \widetilde{\chi}^0$' +
                         '\n'+r'$\tau_{\widetilde{g}}= $'+str(lifetime_days)+
                         ' days\nConstruction time = '+str(construction_time_days)+
                         ' days\nAbsorption time = '+str(2*lifetime_days)+
                         ' days\nDetection time = '+str(30)+
                         ' days\nNumber of RPC hits required = '+str(NRpcsFiring)+
                         ' \nPosition '+str(position)
                )

                plt.tight_layout()
                plt.savefig("plots/sensitivity_"+str(lifetime_days)+'_'+str(construction_time_days)+"_pos"+str(position)+'_nRpcFiring'+str(NRpcsFiring)+'.pdf')
                
def makeLifetimePlots(position):

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    
    for NRpcsFiring in [2]:
        for construction_time_days in [7]:
            for absorption_time_days in [14]:
                for detection_time_days in [30]:
                    for deltam in [3, 300]:

                        xPoints, yPoints = np.meshgrid(np.logspace(0,np.log10(1000)),np.logspace(0,np.log10(1000)))
                        nSigEvents = np.zeros((len(xPoints),len(yPoints)))
                        nBkgEvents = np.zeros_like(nSigEvents)
            
                        for l in range(len(xPoints)):
                            for m in range(len(yPoints)):
                                if xPoints[m,l]!=0.:
                                    nSigEvents[m,l] = nEventsDetectable(
                                        massClass(yPoints[m,l],position),
                                        lifetime_seconds = 3600.*24*xPoints[m,l],
                                        absorption_time_days = absorption_time_days,
                                        construction_time_days = construction_time_days,
                                        detection_time_days = detection_time_days
                                    )

                                    nBkgEvents[m,l] = calcNBkgEvents(
                                        detection_time_days = detection_time_days,
                                        M = NRpcsFiring,
                                        cosmicRate = getMuonsPerSecAboveThresholdBuffered(
                                            deltam/2.,
                                            1. #for 1 m2 since it's already in the function
                                        )
                                    )
                                    #print("for lifetime of "+str(xPoints[m,l])+" days, gluino mass of "+str(yPoints[m,l])+", nSigEvents["+str(m)+","+str(l)+"] is: "+str(nSigEvents[m,l])) 

                        #nSigEvents = np.where(nSigEvents==0, np.nan, nSigEvents)
                        signif = nSigEvents/np.sqrt(nSigEvents+nBkgEvents)
            
                        if deltam == 3:
                            color = 'b'
                            manual_locations = [(200,500),(500,300)]
                            ax.text(8,300,r'$\Delta m = 3$'+' GeV',color=color)
                        elif deltam == 300:
                            color = 'k'
                            manual_locations = [(200,800),(500,500)]
                            ax.text(1.7,700,r'$\Delta m = 300$'+' GeV',color=color)

                        #c = ax.pcolormesh(xPoints, yPoints, signif, cmap='viridis', vmin=0, vmax=5)
                        c = ax.contour(xPoints, yPoints, signif, levels=[3,5], colors=color, linestyles=('dashed','solid'))
                        ax.clabel(c, inline=True, fmt='%1.1f', fontsize='smaller', inline_spacing=17, use_clabeltext=True, manual=manual_locations)
            
                        ax.axis([np.min(xPoints), np.max(xPoints), 10, np.max(yPoints)])
                        ax.set_xlabel(r"$\tau_{\widetilde{g}}$"+" [days]")
                        ax.set_ylabel(r"$m_{\widetilde{g}}$"+" [GeV]")
                        ax.set_xscale('log')
                        ax.set_yscale('log')

                        #fig.colorbar(c, ax=ax,label='$S/\sqrt{S+B}$')

                        #superimpose CMS observed limits (doi:10.1007/JHEP05(2018)127, https://www.hepdata.net/record/ins1645630?version=1&table=Table%2015)
                        #cmsFile = np.genfromtxt('HEPData-ins1645630-v1-Table_15.csv', delimiter=",", names=["x", "y"])
                        #cms = plt.plot(cmsFile['x']/3600/24,cmsFile['y'], color='red')
                        #ax.text(1.2, 700, "CMS 95% CL\nobserved limits", color='red')

    ax.text(2, 20,
            r'$\widetilde{g} \to g \widetilde{\chi}^0$' +
            '\nConstruction time = '+str(construction_time_days)+
            ' days\nAbsorption time = '+str(absorption_time_days)+
            ' days\nDetection time = '+str(30)+
            ' days'+'\nNumber of RPC hits required = '+str(NRpcsFiring)+
            ' \nPosition '+str(position)
    )

    plt.tight_layout()
    plt.savefig("plots/sensitivityVsLifetime_pos"+str(position)+'_absTime'+str(absorption_time_days)+'_nRpcFiring'+str(NRpcsFiring)+'_constrTime'+str(construction_time_days)+ '.pdf')


def printBkgEvents():

    print("number of bkg events is: "+str(int(calcNBkgEvents(detection_time_days = 30)))+" for 30 days detection time and 4 RPCs, 2 of which detected the muon")
    print("number of bkg events is: "+str(int(calcNBkgEvents(M = 3, detection_time_days = 30)))+" for 30 days detection time and 4 RPCs, 3 of which detected the muon")
    print("number of bkg events is: "+str(int(calcNBkgEvents(M = 4, detection_time_days = 30)))+" for 30 days detection time and 4 RPCs, 4 of which detected the muon")
    

#printBkgEvents()
makeplots(0) #comment to use it as a package
#makeplots(1)
makeLifetimePlots(0)
#makeLifetimePlots(1)


