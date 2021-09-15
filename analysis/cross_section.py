import math

import numpy as np

def hasDecayedProb(t, lifetime):
    scaledlt = t*math.log(2)
    return 1. - np.exp( -t/lifetime ) #FIXME: check

class massClass(object):
    def __init__(self, mass : str):
        self.totalefficiency = 0.01
        #in reality check if mass string exists in LUT here
    def eff(self):
        return self.totalefficiency

def totalAbsEfficiency(mass : massClass):
    return mass.eff() #this should be a look-up table based on Jasmine's plots


def lumiInTimeFb(
        t_seconds,
        effective_collision_time = 0.5, #half of the time collisions
        L_int :float = 2.5e34 #on average half of peak lumi
        ):
    
    return t_seconds*effective_collision_time*L_int*1e-39
    

#totalHLXCheck = lumiInTimeFb(3600*24*356*9)
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
    


########### THE STUFF BELOW IS JUST FOR TESTING


def makeplots():
    
    import matplotlib.pyplot as plt
    for lifetime_days in [7., 30., 365]:
        for construction_time_days in [7.]:
        
            x, y = np.meshgrid(np.linspace(lifetime_days/7., 15*lifetime_days, 200), 
                                           np.linspace(lifetime_days/7., 30.*lifetime_days, 200))   
            
            
            z = crossSectionNeededFb(
                massClass('5'),
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
            
            
            plt.savefig("times_"+str(lifetime_days)+'_'+str(construction_time_days)+'.pdf')
    
#makeplots() #comment to use it as a package    