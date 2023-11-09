import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import matplotlib as mpl

mr = [25,250,500,1000,1836]
c = ['b','g','r','c','m','y','k']
l = [(0, (5, 10)),':','-.','--','-']
N_data=250
t = np.linspace(0,6,N_data)

for i in range(len(mr)):

    m=mr[i]

    f = pd.read_csv('ENERGIES_2048_' + str(m)+ '_normalised.csv')

    plt.plot(t,f['Ebn']-f['Ebn'][0],linestyle=l[i],linewidth=0.6,color=c[0],label=r'$\Delta \epsilon^{B}$')
    plt.plot(t,f['Een']-f['Een'][0],linestyle=l[i],linewidth=0.6,color=c[1],label=r'$\Delta \epsilon^{E}$')
    plt.plot(t,f['Eifln']-f['Eifln'][0],linestyle=l[i],linewidth=0.6,color=c[2],label=r'$\Delta \epsilon^{fl}_{i}$')
    plt.plot(t,f['Eefln']-f['Eefln'][0],linestyle=l[i],linewidth=0.6,color=c[3],label=r'$\Delta \epsilon^{fl}_{e}$')
    plt.plot(t,f['Eithn']-f['Eithn'][0],linestyle=l[i],linewidth=0.6,color=c[4],label=r'$\Delta \epsilon^{th}_{i}$')
    plt.plot(t,f['Eethn']-f['Eethn'][0],linestyle=l[i],linewidth=0.6,color=c[5],label=r'$\Delta \epsilon^{th}_{e}$')
    plt.plot(t,f['Etotn']-f['Etotn'][0],linestyle=l[i],linewidth=0.6,color=c[6],label=r'$\Delta \epsilon^{total}$')

#plt.legend()
plt.savefig('all_energies_norm.pdf')
plt.show()






