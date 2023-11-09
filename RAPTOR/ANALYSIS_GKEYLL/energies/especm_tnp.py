import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import raptor
#import sys
#import os
#for dd in ['sims/MR_sims', 'analysis/MR_analysis', 'analysis/init']:
#    if os.path.exists('/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/'):
#        sys.path.insert(0,'/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/')
import lpg
from analysis_engine.statistics.spectra import per_spectra, spectra_base
from analysis_engine.statistics import moments
import pickle


def test_parsevals(ar, fek, L, name=''):
    ce = moments.var(ar, lenn=[L, L])
    fe = np.sum(fek)
    if not np.isclose(ce, fe):
        print(name, ce, fe, ce/fe)
    else:
        print('%s: %s passed!' % (name, ce))





#number of data points
N_data = 250


mr =[25,50,100,250,500,1000,1836] #250,500,1000,1836]
#mr = [25,250,1836,50,500,100,1000]

c = ['b','g','r','c','m','y','k']


l = [(0, (1, 10)),(0, (3, 10, 1, 10)),(0, (5, 10)),':','-.','--','-']



ddd = ['init204825','init204850','init2048100','init2048250','init2048500','init20481000','init20481836']
#ddd = ['init204825','init2048250','init20481836','init204850','init2048500','init2048100','init20481000']

k1, k2, k3 = [],[],[]
fek1, fek2, fek3 = [],[],[]
di = []
g=[]

d={}

r={}

recompute=input('recomput? ')

if recompute == 'y':
    for i in range(len(ddd)):
    
        m = mr[i]
        d_i = np.sqrt(m)
        di.append(d_i)
        L = 8.*np.pi*d_i
        LL = 8. * np.pi
        a = ddd[i]
    
        #loading the specific energies 
        data = raptor.load_data(a, 204, ('bx','by','bz','jex','jey','jez','jix','jiy','jiz','rme','rmi'))
        bx, by, bz = data.bx[:,:,0], data.by[:,:,0], data.bz[:,:,0]             
        bz = bz - np.mean(bz)
        jex, jey, jez = data.jex[:,:,0], data.jey[:,:,0], data.jez[:,:,0]
        jix, jiy, jiz = data.jix[:,:,0], data.jiy[:,:,0], data.jiz[:,:,0]
        rme, rmi = data.rme[:,:,0], data.rmi[:,:,0]
        #print(rme, rmi)
    
        #Calc.i for B field energy spec
        bxs, bys, bzs = (bx*d_i)/(m*0.0023), (by*d_i)/(m*0.0023), (bz*d_i)/(m*0.0023)
        B = (bxs**2 + bys**2 + bzs**2)
    
        #Calc. for ion flow spec
        jexs, jeys, jezs = jex/rme, jey/rme, jez/rme
        jixs, jiys, jizs = jix/rmi, jiy/rmi, jiz/rmi 
        #jexs, jeys, jezs, jixs, jiys, jizs = jex,jey,jez,jix,jiy,jiz
        efl = (jexs**2 + jeys**2 + jezs**2)
        ifl = (jixs**2 + jiys**2 + jizs**2)
    
        k11, fekxB, fekyB, fekzB, fekB = spectra_base.calculate_integrated_spectrum((bxs,bys,bzs), method='periodogram', spec_type='omni', lenn=[L for _ in range(B.ndim)])
    
        test_parsevals(bxs, fekxB, L, 'bx')
        test_parsevals(bys, fekyB, L, 'by')
        test_parsevals(bzs, fekzB, L, 'bz')
    
        fek11 = fekB/L**2 
        k1.append(k11)
        fek1.append(fek11)
        
        k22, fekjex, fekjey, fekjez, fekfle = spectra_base.calculate_integrated_spectrum((jexs,jeys,jezs), method='periodogram', spec_type='omni', lenn=[L for _ in range(efl.ndim)])
    
        test_parsevals(jexs, fekjex, L, 'jex')
        test_parsevals(jeys, fekjey, L, 'jey')
        #test_parsevals(jezs, fekjez, L, 'jez')
    
        ## NOTE: The different normalization
        fek22 = fekfle/L**2
        k2.append(k22)
        fek2.append(fek22)
    
        k33, fekjix, fekjiy, fekjiz, fekfli = spectra_base.calculate_integrated_spectrum((jixs,jiys,jizs), method='periodogram', spec_type='omni', lenn=[L for _ in range(ifl.ndim)])
    
        test_parsevals(jixs, fekjix, L, 'jix')
        test_parsevals(jiys, fekjiy, L, 'jiy')
        #test_parsevals(jizs, fekjiz, L, 'jiz')
    
        fek33 = fekfli/L**2
        k3.append(k33)
        fek3.append(fek33)
        
        #d[m] = {'bx':bxs,'by':bys,'bz':bzs,'jex':jex,'jey':jey,'jez':jez,'rme':rme,'jix':jix, 'jiy':jiy,'jiz':jiz,'rmi':rmi} 
        d[i] = {'k':k11,'fek1':fek11,'fek2':fek22,'fek3':fek33}
    pickle.dump(d,open('omni-spectra.pkl','wb'))
else:
    d = pickle.load(open('omni-spectra.pkl','rb'))
    
print(d)
print(len(d[0]['fek1']))
print(len(d[0]['k']))

#np.save('omni-spectra-data.npy',d)

plt.rcParams.update({'font.size': 9})
#plt.rc('text', usetex=True)

fig, ax = plt.subplots(3, 1, sharex=True,figsize=(3.5,6))
fig.subplots_adjust(hspace=0.0)#5)


for h in range(len(ddd)):
    k1=k2=k3=d[h]['k']
    fek1a=d[h]['fek1']
    fek2a=d[h]['fek2']
    fek3a=d[h]['fek3']
    
    print(len(d[h]['fek1']))
    print(len(fek1a))
    print(len(k1))

    fm1 = np.max(fek1a)
    fm2 = np.max(fek2a)
    fm3 = np.max(fek3a)
    m = mr[h]
    d_i = np.sqrt(m)

    print(len(fek1a/fm1))
    print(len(k1*d_i))

    ax[0].loglog(k1*d_i,fek1a/fm1,linestyle=l[h],linewidth='0.6',color=c[h],label=str(m))
    ax[0].axvline(d_i, color=c[h], linestyle=':',linewidth=0.4,alpha=0.8)
    ax[1].loglog(k2*d_i,fek2a/fm2,linestyle=l[h],linewidth='0.6',color=c[h])
    ax[1].axvline(d_i, color=c[h], linestyle=':',linewidth=0.4,alpha=0.8)
    ax[2].loglog(k3*d_i,fek3a/fm3,linestyle=l[h], linewidth='0.6',color=c[h])
    ax[2].axvline(d_i, color=c[h], linestyle=':',linewidth=0.4,alpha=0.8)

#labels
ax[0].text(0.3,0.00015,r'$\epsilon^{B}(k)$')
#ax[0].set_title(r'- 25  -.- 50  -- 100')
#ax[0].set_title(r'Omni-Spectra at $\tau=4.9\tau_{0}$')

ax[1].text(0.3,0.0015,r'$\epsilon^{ve}(k)$')

ax[2].text(0.3,0.00015,r'$\epsilon^{vi}(k)$')
ax[2].set_xlabel(r'$k d_i$')

# Put a legend below current axis
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=4)

#kdi=1 lines
ax[0].axvline(1, color='gray', linestyle=':')
ax[1].axvline(1, color='gray', linestyle=':')
ax[2].axvline(1, color='gray', linestyle=':')

#ticks inward
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




plt.savefig('enery_spectra_multi_late_times.pdf',bbox_inches='tight',pad_inches=0.05)

plt.show()
























