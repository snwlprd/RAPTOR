import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import raptor

mr = [25,50,100,250,500,1000,1836]
dd = ['init204825','init204850','init2048100','init2048250','init2048500','init20481000','init20481836']
N_data = 250

t = np.linspace(0,6,250)

for i in range(len(mr)):
    print(i)

    m = mr[i]
    d=dd[i]
    #loading Energy data 
    f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_'+str(m)+'.csv')

    #loading PiD, pTh, pgu
    g = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/pgu/PI_D_VALUES_'+str(m)+'.csv')
    
    #defining data to plot
    PiD_ions, PiD_elec = g['PiD_ions'], g['PiD_elec']
    pTh_ions, pTh_elec = g['pTh_ions'], g['pTh_elec'] 
    pgu_ions, pgu_elec = g['pgu_ions'], g['pgu_elec']

    norm = g['norm']
    #normm = f['norm']

    #j dot E
    jde, sjde = [], []

    for j in range(0,N_data):
        print(j)

        l_jde = raptor.j_dot_E(d,m,j)
        a_jde = np.mean((l_jde))
        jde.append(a_jde)

        if j ==0:
            sjde.append(jde[j])
        else: 
            sum_jde = jde[j] + sjde[j-1]
            sjde.append(sum_jde)

    sjde_np = np.array(sjde)
    
    #calculating tau_naught
    tau_zero = (8*np.pi*np.sqrt(m))/(2*np.pi*0.0023)

    #calculating dt
    dt = 6*tau_zero/250

    c_sjde_np = sjde_np*dt  

    #plotting 
    plt.plot(t,(PiD_ions)/(norm),label=r'$\Pi D$')
    plt.plot(t,(pTh_ions)/(norm),label=r'$p \Theta$')
    plt.plot(t,(pgu_ions)/(norm),label=r'$(P \cdot \nabla) \cdot u$')
    plt.plot(t,(f['Eith']-f['Eith'][0])/(norm),linestyle='--',label=r'$\Delta \epsilon^{th}_{i}$')
    plt.legend()
    plt.title(r'Ions with $m_{i}/m_{e}$ = ' + str(m))
    plt.xlabel(r'Eddy turnover time, $\tau_{0}$')
    plt.ylabel(r'Cumulative $\Pi D$, Cumulative $p \Theta$')
    plt.savefig('cumulative_PiD_pTh_vs_Eith_'+str(m)+'.pdf')
    #plt.show()
    plt.clf()
    
    plt.plot(t,(PiD_elec)/(norm),label=r'$\Pi D$')
    plt.plot(t,(pTh_elec)/(norm),label=r'$p \Theta$')
    plt.plot(t,(pgu_elec)/(norm),label=r'$(P \cdot \nabla) \cdot u$')
    #plt.plot(t,(c_sjde_np)/(norm),linestyle='-.',label=r'$j \cdot E$')
    plt.plot(t,(f['Eeth']-f['Eeth'][0])/(norm),linestyle='--',label=r'$\Delta \epsilon^{th}_{e}$')
    plt.legend()
    plt.title(r'Electrons with $m_{i}/m_{e}$ = ' + str(m))
    plt.xlabel(r'Eddy turnover time, $\tau_{0}$')
    plt.ylabel(r'Cumulative $\Pi D$, Cumulative $p \Theta$')
    plt.savefig('cumulative_PiD_pTh_vs_Eeth_'+str(m)+'.pdf')
    #plt.show()
    plt.clf()


