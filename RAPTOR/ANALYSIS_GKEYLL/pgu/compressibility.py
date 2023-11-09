import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import raptor

mr = [25,50,100,250,500,1000,1836]

N_data = 250

t = np.linspace(0,6,250)

#r_e, r_i = [], []

plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots()

for i in range(len(mr)):
    print(i)

    m = mr[i]
    
    #loading Energy data 
    #f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_test/analysis/energies/ENERGIES_TEST_'+str(m)+'.csv')

    #loading PiD, pTh, pgu
    g = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_test/analysis/pid/PI_D_VALUES_'+str(m)+'.csv')
    
    #defining data to plot
    PiD_ions, PiD_elec = g['PiD_ions'], g['PiD_elec']
    pTh_ions, pTh_elec = g['pTh_ions'], g['pTh_elec'] 
    #pgu_ions, pgu_elec = g['pgu_ions'], g['pgu_elec']

    norm = g['norm']
    
    #calculating tau_naught
    tau_zero = (8*np.pi*np.sqrt(m))/(2*np.pi*0.0023*0.1)

    #calculating dt
    dt = 6*tau_zero/250

    c_i = np.mean(pTh_ions[167:250]/PiD_ions[167:250]) #np.mean(pTh_ions[167:250])/np.mean(PiD_ions[167:250])
    c_e = np.mean(pTh_elec[167:250]/PiD_elec[167:250]) #np.mean(pTh_elec[167:250])/np.mean(PiD_elec[167:250]) 
    
    ax.loglog(m,c_i,'+',color='k')
    ax.loglog(m,c_e,'x',color='k')

#h = linspace()

ax.tick_params(axis='y',direction='in')
ax.tick_params(axis='x',direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(which='minor', direction='in')

ax.text(25,1.4,'x-electrons',color='k')
ax.text(25,1.0,'+-ions',color='k')

#plotting 
plt.plot()
#plt.legend()
#plt.title(r'Compressibility for $4\tau_{0} \leq \tau \leq 6\tau_{0}$')
plt.xlabel(r'$m_{i}/m_{e}$')
plt.ylabel(r'$<p\theta$ / $\Pi D>$')
plt.savefig('compresibility.pdf',bbox_inches='tight',pad_inches=0.05)
plt.show()

