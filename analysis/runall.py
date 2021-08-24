from multiprocessing import Pool
import os
import time

parttemp='''      1000021   {mass}.0     # ~g
      1000022   {massn}     # Neutralino
      1000993   {mass}.700   # ~g_glueball
      1009213   {mass}.700   # ~g_rho+
      1009113   {mass}.700   # ~g_rho0
      1093122   {mass}.700   # ~g_Lambda0
     -1009213   {mass}.700   # ~g_rho-
'''


def worker(no):
        
    beta,fileno = betas[no],no
    os.system(" ./exampleB4a -b {beta} -o {offset} -f {file} -m  {macro} -a {add}".format(
        beta=beta,
        file=fileno,
        offset=fileno,
        macro="run2.mac",
        add = add
        ))
    return True
    
def run_mass(mass: int):
    
    global betas
    global add
    
    with open("particles.txt","w") as f:
        f.write(parttemp.format(mass=mass,massn=mass/2.))

    betas = range(1*1,1*99)
    betas = [b/100. for b in betas]
    add = "_m"+str(mass)
    
    
    
    #maxpar = min(len(betas),30)
    print('running',len(betas),'jobs')
    p = Pool()
    ret = p.map(worker, range(len(betas)))
    
    

for mass in [5,10,20,40,80,160,320,640,1280,2560]:
    run_mass(mass) 
    time.sleep(1) #some time before next mass point 
    



