import numpy as np
import raptor
import matplotlib.pyplot as plt
import pandas as pd

N = 250
mr = [25,50,100,250,500,1000,1836] #,250,500,1000,1836]
jz = np.zeros((len(mr),N))
c = ['b','g','r','c','m','y','k']
l = [(0, (1, 10)),(0, (3, 10, 1, 10)),(0, (5, 10)),':','-.','--','-']
d = ['init204825','init204850','init2048100','init2048250','init2048500','init20481000','init20481836']

plt.rcParams.update({'font.size': 9})


t=np.linspace(0,6,250)

'''
for i in range(len(mr)):
    print(mr[i])
    m = mr[i]
    a = d[i]
    for j in range(0,N):
        print(j)
        jxl, jyl, jzl = raptor.j_xyz(a,m,j)
        jz[i][j] = np.sqrt(np.mean(jzl**2))

    #make pandas data frame and save energies
    Eu  = {'t':t, 'rms':jz[i]}
    uf = pd.DataFrame(Eu)

    uf.to_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/jz/RMS_JZ_{}.csv'.format(str(m)))
'''
fig,ax=plt.subplots()

for i in range(len(mr)):
    print(mr[i])
    m=mr[i]

    f = pd.read_csv('RMS_JZ_' +str(m)+ '.csv')
    jz[i]=f['rms']

ax.tick_params(axis='y',direction='in')
ax.tick_params(axis='x',direction='in')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(which='minor', direction='in')


mxm = np.max(jz[6])

for k in range(len(mr)):
    ax.plot(t,jz[k]/mxm,linewidth=0.6,color=c[k],linestyle=l[k],label=str(mr[k]))
    



plt.legend()
plt.xlabel(r'$\tau_{0}$')
plt.ylabel(r'rms( $j_{z}$ )')
plt.savefig('rms-jz.pdf',bbox_inches='tight',pad_inches=0.05)
plt.show()
