import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import raptor

#number of data points
N_data = 250

#defining the Mass Ratio 
mr = [25,50,100,250,500,1000,1836]
c = ['b','g','r','c','m','y','k']

aa, bb, aaa, bbb =[], [], [], []


for i in range(len(mr)):
    print(i)
    
    m = mr[i]

    #loading in data from chosen file
    f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_' + str(m) + '.csv')

    #Average derivative value for all the energies
    a=f.diff(periods=1,axis=0)
    qi = np.mean(a['Eith'][166:])
    qe = np.mean(a['Eeth'][166:])
    qi_qe = qi/qe
    qpq = qi + qe
    print(qi)
    print(qe)
    aaa.append(qi)
    bbb.append(qe)
    aa.append(qi_qe)
    bb.append(qpq)
'''
for j in range(len(mr)):
    m = mr[j]
    plt.plot(bb[j],aa[j],'x',color=c[j])
    plt.ylabel(r'$<Q_{i}>/<Q_{e}>$')
    plt.xlabel(r'$m_{i}/m_{e}$')
    plt.title(r'$<Q_{i}>/<Q_{e}>$ vs $m_{i}/m_{e}$')
    plt.savefig('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_sims/comparison_analysis/qi_qe/turb_region/qi_qe_vs_qpq.svg')
plt.clf()
'''
#font and font size
plt.rcParams.update({'font.size': 9})
#plt.rc('text', usetex=True)

fig, ax = plt.subplots(4,1,sharex=True,figsize=(3.5,6))
fig.subplots_adjust(hspace=0.1)

for k in range(len(mr)):
    m = mr[k]
    ax[0].semilogx(m,bbb[k]/bbb[6],'x',linestyle=':',color='k')
    #ax.set_xscale('log')
    #ax[0].set_ylabel(r'$<Q_{i}>$')
    #ax[0].set_title(r'$<Q_{i}/<Q_{e}>$, $<Q_{i}>$ and $<Q_{e}>$ vs $m_{i}/m_{e}$ for $4 \leq \tau_{0} \leq 6$')
    ax[1].semilogx(m,aaa[k]/aaa[6],'x',color='k') #c[2])
    #ax[1].set_ylabel(r'$<Q_{e}>$')
    ax[2].semilogx(m,bb[k]/bb[6],'x',color='k')#c[3])
    #ax[2].set_ylabel(r'$<Q_{i}>/<Q_{e}>$')
    ax[3].semilogx(m,aa[k]/aa[6],'x',color='k')#c[4])
    #ax[3].set_ylabel(r'$<Q_{i}> + <Q_{e}>$')
    ax[3].set_xlabel(r'$m_{i}/m_{e}$')
    
    #ax[0].title(r'$<Q_{i}> and <Q_{e}>$ vs $m_{i}/m_{e}$')

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
ax[0].text(25,0.8,r'$<Q_{e}>$')
ax[1].text(25,0.8,r'$<Q_{i}>$')
ax[2].text(25,0.8,r'$<Q_{i}> + <Q_{e}>$')
ax[3].text(25,0.95,r'$<Q_{i}>/<Q_{e}>$')

# Put a legend above current axis
#ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=4)





plt.show()
fig.savefig('qi_qe_combination_comp_v1.pdf')
plt.clf()

'''
for h in range(len(mr)):
    m = mr[h]
    plt.plot(m,aa[h],'x',color=c[h])
    plt.plot(m,bb[h],'+',color=c[h])
    plt.ylabel(r'$<Q_{i}>/<Q_{e}>$, $<Q_{i}>+<Q_{e}>$')
    plt.xlabel(r'$m_{i}/m_{e}$')
    plt.title(r'$<Q_{i}>/<Q_{e}>$(x), $<Q_{i}>+<Q_{e}>$(+) vs $m_{i}/m_{e}$')
    plt.savefig('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_sims/comparison_analysis/qi_qe/turb_region/qi_qe_qpq_mr.svg')
plt.clf()
'''

















