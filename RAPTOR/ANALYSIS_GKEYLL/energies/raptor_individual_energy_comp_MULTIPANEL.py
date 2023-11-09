import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import matplotlib as mpl

mr = [25,50,100,250,500,1000,1836]
c = ['b','g','r','c','m','y','k']
l = [(0, (1, 10)),(0, (3, 10, 1, 10)),(0, (5, 10)),':','-.','--','-']
N_data=250
t = np.linspace(0,6,N_data)
ef=[]

plt.rcParams.update({'font.size': 9})
#plt.rc('text', usetex=True)

fig, ax = plt.subplots(4, 1, sharex=True,figsize=(3.5,6))
fig.subplots_adjust(hspace=0.1)

for i in range(len(mr)):

    m=mr[i]

    f = pd.read_csv('ENERGIES_2048_' + str(m)+ '_normalised.csv')

    ax[0].plot(t,f['Ebn']-f['Ebn'][0],linestyle=l[i],linewidth=0.6,color=c[i],label=str(m))
    #ax[0].set_ylabel(r'$\Delta \epsilon^{B}$')
    #ax[0].set_title(r'$\Delta \epsilon$ for different $m_{i}/m_{e}$')
    ax[1].plot(t,f['Eifln']-f['Eifln'][0],linestyle=l[i],linewidth=0.6,color=c[i])
    #ax[1].set_ylabel(r'$\Delta \epsilon^{fl}_{i}$')
    ax[2].plot(t,f['Eithn']-f['Eithn'][0],linestyle=l[i],linewidth=0.6,color=c[i])
    #ax[2].set_ylabel(r'$\Delta \epsilon^{th}_{i}$')
    #ax[2].set_xlabel(r'$Eddy Turnover time, \tau_{0}$')
    ax[3].plot(t,f['Eethn']-f['Eethn'][0],linestyle=l[i],linewidth=0.6,color=c[i])
    #ax[3].set_ylabel(r'$\Delta \epsilon^{th}_{e}$')
    ax[3].set_xlabel(r'$\tau_{0}$')

#set ticks inward
ax[0].tick_params(axis='y',direction='in')
ax[0].tick_params(axis='x',direction='in')
ax[0].xaxis.set_ticks_position('both')
ax[0].yaxis.set_ticks_position('both')
ax[0].tick_params(which='minor', direction='in')

ax[1].tick_params(axis='y',direction='in')
ax[1].tick_params(axis='x',direction='in')
ax[1].xaxis.set_ticks_position('both')
ax[1].yaxis.set_ticks_position('both')
ax[1].tick_params(which='minor', direction='in')

ax[2].tick_params(axis='y',direction='in')
ax[2].tick_params(axis='x',direction='in')
ax[2].xaxis.set_ticks_position('both')
ax[2].yaxis.set_ticks_position('both')
ax[2].tick_params(which='minor', direction='in')

ax[3].tick_params(axis='y',direction='in')
ax[3].tick_params(axis='x',direction='in')
ax[3].xaxis.set_ticks_position('both')
ax[3].yaxis.set_ticks_position('both')
ax[3].tick_params(which='minor', direction='in')

#set label inside figure
ax[0].text(0.2,0.125,r'$\Delta \epsilon^{B}$')
ax[1].text(0.2,-0.125,r'$\Delta \epsilon^{fl}_{i}$')
ax[2].text(0.2,0.225,r'$\Delta \epsilon^{th}_{i}$')
ax[3].text(0.2,0.07,r'$\Delta \epsilon^{th}_{e}$')

# Put a legend above current axis
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=4)

plt.savefig('energies_comp_norm_2048.pdf',bbox_inches='tight',pad_inches=0.05)
plt.show()





'''
Reverse list
colors=[][::-1]
colors=['b','k',...][::-1]
'''







