import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
for dd in ['WorkSpace', 'WorkSpace/python', 'AJGAR/TurbAn']:
    if os.path.exists(os.environ['HOME']+'/'+dd):
        sys.path.insert(0,os.environ['HOME']+'/'+dd+'/')
for dd in ['sims/MR_2048','sims/MR_sims', 'analysis/MR_analysis', 'analysis/init']:
    if os.path.exists('/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/'):   
        sys.path.insert(0,'/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/')
import lpg


 
#number of data points
N_data          = 250

r = {}

Eb   = np.empty(N_data)
Ee   = np.empty(N_data)
Eeth = np.empty(N_data)
Eith = np.empty(N_data)
Eefl = np.empty(N_data)
Eifl = np.empty(N_data)
Etot = np.empty(N_data)

#grid size
grid = [50,100] #[25,250,500,1000,1836]

#calculating and saving energies
for i in range(len(grid)):
    print("mr{}".format(grid[i]))

    r[grid[i]]=lpg.lpg('init2048'+str(grid[i]))
    r[grid[i]].vars2load(['bx', 'by', 'bz', 'ex', 'ey', 'ez', 'jex', 'jey', 'jez', 'jix', 'jiy', 'jiz', 'pexx', 'peyy', 'pezz', 'pixx', 'piyy', 'pizz', 'rme', 'rmi'])

    for j in range(0, N_data):
        print("time step {}".format(j))
        r[grid[i]].loadslice(j)

        #calculating field energies
        bx, by, bz = r[grid[i]].bx, r[grid[i]].by, r[grid[i]].bz
        dbx, dby, dbz = bx - np.mean(bx), by - np.mean(by), bz - np.mean(bz)
        B       = dbx**2 + dby**2 + dbz**2
        
        E       = r[grid[i]].ex**2 + r[grid[i]].ey**2 + r[grid[i]].ez**2
        Eb[j]   = 0.5*np.mean(B)     
        Ee[j]   = 0.5*np.mean(E)    

        #calculating flow energies
        efl     =  (r[grid[i]].jex**2 + r[grid[i]].jey**2 + r[grid[i]].jez**2)/r[grid[i]].rme
        ifl     =  (r[grid[i]].jix**2 + r[grid[i]].jiy**2 + r[grid[i]].jiz**2)/r[grid[i]].rmi
        Eefl[j] = 0.5*np.mean(efl)     
        Eifl[j] = 0.5*np.mean(ifl)     
 
        #calculating thermal energies
        eth     = r[grid[i]].pexx + r[grid[i]].peyy + r[grid[i]].pezz
        ith     = r[grid[i]].pixx + r[grid[i]].piyy + r[grid[i]].pizz
        Eeth[j] = (0.5*np.mean(eth)) - Eefl[j]    #2*Eefl[j]
        Eith[j] = (0.5*np.mean(ith)) - Eifl[j]    #2*Eifl[j]

        #calculating total energy
        #Etot[j] = Eb[j] + Ee[j] + Eefl[j] + Eifl[j] + Eeth[j] + Eith[j]
        Etote = 0.5*eth     
        Etoti = 0.5*ith
        Etot[j] = np.mean(Etote) + np.mean(Etoti) + Eb[j] + Ee[j]

    #timestep
    t = np.linspace(0,6,250)

    #make pandas data frame and save energies 
    Eu  = {'t':t, 'Eb':Eb, 'Ee':Ee, 'Eefl':Eefl, 'Eifl':Eifl, 'Eeth':Eeth, 'Eith':Eith, 'Etot':Etot}
    uf = pd.DataFrame(Eu)

    uf.to_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_{}.csv'.format(str(grid[i])))

    #finding energy in dB
    r[grid[i]].loadslice(0)
    bzero =np.sqrt(r[grid[i]].bx**2 + r[grid[i]].by**2 + r[grid[i]].bz**2)
    
    r[grid[i]].loadslice(1)
    bone =np.sqrt(r[grid[i]].bx**2 + r[grid[i]].by**2 + r[grid[i]].bz**2)

    EdB = Eb[1] + Eb[0] - 0.5*np.mean(2*bzero*bone)

    #E total 0
    En = Eb[0] + Eefl[0] + Eifl[0]
    
    #normalising 
    Ebn = Eb/En
    Een = Ee/En
    Eefln = Eefl/En
    Eifln = Eifl/En
    Eethn = Eeth/En
    Eithn = Eith/En
    Etotn = Etot/En
   
    #make pandas data frame and save normalised energies 
    E  = {'t':t, 'Ebn':Ebn, 'Een':Een, 'Eefln':Eefln, 'Eifln':Eifln, 'Eethn':Eethn, 'Eithn':Eithn, 'Etotn':Etotn}
    df = pd.DataFrame(E)
    
    df.to_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_{}_normalised.csv'.format(str(grid[i])))




