
import numpy as np

'''
Just a small module to get expected background events from a muon spectrum 
convolved with the maximal deposit in the detector volume
'''

def getLArEnergy(muon_energy):
    return 0.01*muon_energy #to be refined, but seems to be accurate < 10TeV
def getMuonEnergy(LAr_energy):
    return 100.*LAr_energy

def MuonSpectrum(p):
    '''
    vertical muon spectrum - horizontal will just be less energy so should be good upper bound
    from https://user-web.icecube.wisc.edu/~aya/simulation/prompt/ref/review/hep-ph9803488-BMNSST.pdf
    Table 2
    
    per (cm2 * s * sr)
    '''
    
    p = np.array(p)
    p = np.expand_dims(p,axis=1)
    
    #calculate all
    C  = np.array([[2.95e-3, 1.781e-2, 1.435e1, 1e3]] )
    g0 = np.array([[0.3061, 1.7910, 3.6720, 4.]     ] )
    g1 = np.array([[1.2743, 0.3040, 0., 0.]         ] )
    g2 = np.array([[-0.2630, 0., 0., 0.]            ] )
    g3 = np.array([[0.0252, 0., 0., 0.]             ] )
    
    rangeslo = np.array([[1., 9.2765e2, 1.5878e3, 4.1625e5]    ])
    rangeshi = np.array([[9.2765e2, 1.5878e3, 4.1625e5, 1e18]  ])
    #np.where p.. then pick
    #print(p.shape,g0.shape)
    
    D = C * p **(-( g0 + g1*np.log10(p) + g2*np.log10(p)**2 + g3*np.log10(p)**3 ))
    
    selector = np.logical_and(p>=rangeslo, p<rangeshi)
    D = D[selector]
    return D
    


def getMuonsPerSecAboveThreshold(
        min_detection_energy : float,
        detector_area_in_sqm : float,
                      ):
    '''
    numerically integrate over spectrum and energy 
    (because I am lazy)
    
    Includes 2pi steradian, assuming we only get muons form the top, but in turn assuming the flux is
    as energetic from the top as all other sides - so this is still conservative
    '''
    min_muon_energy = getMuonEnergy(min_detection_energy)
    
    def sum_from_to(fromenergy, toenergy,resolution):
        p = np.arange(fromenergy,toenergy,step=resolution)
        N = MuonSpectrum(p)
        return np.sum(N*resolution)
    
    #in steps
    all = 0
    if min_muon_energy < 50:
        all += sum_from_to(min_muon_energy, 50., 0.001)
    if min_muon_energy < 250:
        all += sum_from_to(max(50.,min_muon_energy), 250., 0.1)
    if min_muon_energy < 1000:
        all += sum_from_to(max(250.,min_muon_energy), 1000., 1.)
        
    all += sum_from_to(max(1000.,min_muon_energy), 1e7, 10.)
    N=all          
    sqcm = 1e4*detector_area_in_sqm
    return 2.*3.1415*sqcm*N
    

muonlut={}
def getMuonsPerSecAboveThresholdBuffered(
        min_detection_energy : float,
        detector_area_in_sqm : float = 1.,
        ):
    global muonlut
    if min_detection_energy in muonlut.keys():
        return muonlut[min_detection_energy]*detector_area_in_sqm
    else:
        n = getMuonsPerSecAboveThreshold(min_detection_energy,1.)
        muonlut[min_detection_energy] = n
        return n*detector_area_in_sqm
    


def test_spectrum():
    
    import matplotlib.pyplot as plt
    
    p = 10.**(np.arange(1,5,step=0.01))
    
    plt.plot(p, p**3*MuonSpectrum(p))#reproduces the spectrum plot from paper
    plt.yscale('log')
    plt.xscale('log')
    plt.show()



def loadLUT():
    import pickle
    global muonlut
    with open('cosmics_lut.pkl','rb') as f:
        muonlut = pickle.load(f)
        
def createLUT():
    import tqdm
    import pickle
    
    xPoints, yPoints = np.meshgrid(np.logspace(0,np.log10(3000),num=400),np.logspace(0,np.log10(3000),num=400))
    for j in tqdm.tqdm(range(len(xPoints))):
        for i in range(j+1):
            deltam = xPoints[i,j] - yPoints[i,j]
            getMuonsPerSecAboveThresholdBuffered(
                                    deltam/2.,
                                    1. #for 1 m2 since it's already in the function
                                    )
    with open('cosmics_lut.pkl','wb') as f:
        pickle.dump(muonlut,f)
 
#createLUT()  
loadLUT()
#print(getMuonsPerSecAboveThreshold(3, 16))








