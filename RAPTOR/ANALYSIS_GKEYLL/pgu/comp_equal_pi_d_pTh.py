import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import raptor

mr = [25,50,100,250,500,1000,1836]
dd = ['init204825','init204850','init2048100','init2048250','init2048500','init20481000','init20481836']
N_data = 250

t = np.linspace(0,6,250)

plt.rcParams.update({'font.size': 9})
plt.rcParams["figure.figsize"] = (3.5,1.75)

ccr=np.zeros(7)
print(ccr)
for i in range(len(mr)):
    print(i)

    m = mr[i]
    d=dd[i]
    #loading Energy data 
    #f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_'+str(m)+'.csv')

    #loading PiD, pTh, pgu
    g = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/pgu/PI_D_VALUES_'+str(m)+'.csv')
    
    #defining data to plot
    PiD_elec = g['PiD_elec']
    pTh_elec = g['pTh_elec'] 
    

    norm = g['norm']
    #normm = f['norm']

    #find cross over point and index
    diff = abs(PiD_elec[40:] - pTh_elec[40:])
    cross = np.where(diff == np.amin(diff))[0]
    cr = cross[-1]
    print(cr)

    #calculating tau_naught
    tau_zero = (8*np.pi*np.sqrt(m))/(2*np.pi*0.0023)

    #calculating dt
    dt = 6*tau_zero/250

    ccr[i] = ((cr+40)/250)*6

print(len(ccr))
#plotting 
plt.plot(mr,ccr,marker='x',linestyle=':',color='k')
plt.xscale('log')
plt.ylim((0,6))
    #plt.plot(t,PiD_elec - pTh_elec)
    #plt.legend()
plt.xlabel(r'$m_i/m_e$')
plt.ylabel(r'$\tau_0$')
    #plt.savefig('equal_PiD_pTh_vs_Eith_'+str(m)+'.pdf')
    #plt.show()
    #plt.clf()

#plt.show()
plt.savefig('equal_PiD_pTh_vs_Eith.pdf',bbox_inches='tight',pad_inches=0.05)
plt.show()
