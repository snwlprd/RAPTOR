import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import matplotlib as mpl

mr = [25,250,500,1000,1836]
c = ['b','g','r','c','m','y','k']
l = [(0, (5, 10)),':','-.','--','-'] #[(0, (1, 10)),(0, (3, 10, 1, 10)),(0, (5, 10)),':','-.','--','-']
N_data=250
t = np.linspace(0,6,N_data)
sim = ['ENERGIES_2048_25']

for i in range(len(mr)):

    m=mr[i]
    
    f = pd.read_csv('ENERGIES_2048_'+str(mr[i])+ '_normalised.csv')

    thm = (f['Eithn']-f['Eithn'][0])+(f['Eethn']-f['Eethn'][0])
    flw = (f['Eifln']-f['Eifln'][0])+(f['Eefln']-f['Eefln'][0])+(f['Ebn']-f['Ebn'][0])
    athm, aflw = np.mean(thm[167:250]), np.mean(flw[167:250])

    plt.plot(t,f['Ebn']-f['Ebn'][0],linestyle=l[i],linewidth=0.6,color=c[0])#,label=r'$\Delta \epsilon^{B}$')
    plt.plot(t,f['Een']-f['Een'][0],linestyle=l[i],linewidth=0.6,color=c[1])#,label=r'$\Delta \epsilon^{E}$')
    plt.plot(t,f['Eifln']-f['Eifln'][0],linestyle=l[i],linewidth=0.6,color=c[2])#,label=r'$\Delta \epsilon^{fl}_{i}$')
    plt.plot(t,f['Eefln']-f['Eefln'][0],linestyle=l[i],linewidth=0.6,color=c[3])#,label=r'$\Delta \epsilon^{fl}_{e}$')
    plt.plot(t,f['Eithn']-f['Eithn'][0],linestyle=l[i],linewidth=0.6,color=c[4])#,label=r'$\Delta \epsilon^{th}_{i}$')
    plt.plot(t,f['Eethn']-f['Eethn'][0],linestyle=l[i],linewidth=0.6,color=c[5])#,label=r'$\Delta \epsilon^{th}_{e}$')
    plt.plot(t,f['Etotn']-f['Etotn'][0],linestyle=l[i],linewidth=0.6,color=c[6],label= str(m))
    #plt.plot(t,thm,linestyle=l[i],linewidth=0.6,color=c[0])
    #plt.plot(t,flw,linestyle=l[i],linewidth=0.6,color=c[1])
    plt.text(0.9,0.8,str(m) + ':' )
    #plt.semilogx(m,athm,'x',color=c[0])
    #plt.semilogx(m,abs(aflw),'.',color=c[1])
   
#plt.text(0.5,0.2,r'$\Delta \epsilon^{thermal}$') 
#plt.text(0.5,-0.2,r'$\Delta \epsilon^{fluctuations}$')
plt.text(0,0.3,'$\Delta \epsilon^{fl}_{i}$:red, $\Delta \epsilon^{fl}_{e}$:cyan')
plt.text(0,0.25,'$\Delta \epsilon^{th}_{i}$:purple, $\Delta \epsilon^{th}_{e}$:yellow')
plt.text(0,0.2,'$\Delta \epsilon^{tot}$:black')
plt.text(0,-0.25,'. . .:25, - . -:50') 
plt.text(0,-0.3,'- - -:100, ......:250')
plt.text(0,-0.35,' -.-.-.:500, ------:1000:  - :1836')
plt.ylabel(r'Normalised Energy Change, $\Delta \epsilon / (\Delta \epsilon^{B} + \Delta \epsilon^{fl}_{i} + \Delta \epsilon^{fl}_{e})$')
plt.xlabel(r'Eddy Turnover Time, $\tau_{0}$')
#plt.xlabel(r'$m_{i} / m_{e}$')
#plt.title(r'average energy for $4\tau_{0} \leq \tau \leq 6\tau_{0}$')
#plt.legend()
plt.savefig('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_test_workload/analysis/energies/energies_norm.pdf')
plt.show()






