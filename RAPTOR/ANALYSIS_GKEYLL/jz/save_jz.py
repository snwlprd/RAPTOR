import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import raptor
import pandas as pd

mr=[1000,1836] #[25,50,100,250,500,1000,1836]
N_data=250
lx, ly, lz = 2048, 2048, 1
#jz = np.zeros((len(mr),lx,ly,lz))
#jx = np.zeros((len(mr),lx,ly,lz))
#jy = np.zeros((len(mr),lx,ly,lz))
L = 8*np.pi 
ex, ey = L/2, L/2
#krts=np.zeros(len(mr))
dd = ['init20481000','init20481836']#['init204825','init204850','init2048100','init2048250','init2048500','init20481000','init20481836']
#abcd = ['a','b','c','d']
#jz = np.empty(N_data)
#jx = np.empty(N_data)
#jy = np.empty(N_data)

t = np.linspace(0,6,250)

for i in range(len(mr)):
    print(i)

    d = dd[i]
    m = mr[i]
    jx, jy, jz = [], [], []

    for j in range(0,N_data):
        print(j)
        jxl, jyl, jzl = raptor.j_xyz(d,m,j)

        jz.append(jzl)
        jx.append(jxl)
        jy.append(jyl)

    #make pandas data frame and save energies
    Eu  = {'t':t, 'jx':jx, 'jy':jy, 'jz':jz}
    uf = pd.DataFrame(Eu)

    uf.to_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/jz/JZ_2048_{}.csv'.format(str(m)))


    #find kurtosis
    #krts[i] = np.mean(jz[i]**4)/(np.mean(jz[i]**2)**2)

#print(np.max(jz))
'''
minmax=[np.min(jz),np.max(jz)]

fig,ax=plt.subplots(1,4,sharey='row',figsize=(10,2))

alpha = np.max(jz)*1.5
print(alpha)

for j in range(len(mr)):
    m=mr[j]
    pilaat=ax[j].imshow(jz[j]/alpha,cmap='seismic',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    ax[j].text(-0.9*ex,0.8*ey,str(m))
    #ax[j].text(0.2,-0.78*ey,'('+str(abcd[j])+')')
    ax[j].set_xlabel('x $(d_{i})$')
    ax[j].add_artist(lines.Line2D([0.406*ex,0.406*ex],[-ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[j].add_artist(lines.Line2D([-0.102*ex,-0.102*ex],[-ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[j].add_artist(lines.Line2D([-0.102*ex,0.406*ex],[-ey,-ey],linestyle=':',linewidth=1,color='black'))
    ax[j].add_artist(lines.Line2D([-0.102*ex,0.406*ex],[-0.829*ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))

ax[0].set_ylabel('y $(d_{i})$')
fig.colorbar(pilaat,ax=ax[:],pad=0.02)
#fig.suptitle(r'$j_{z}/(2.28\times10^{-2})$ at $\tau=4.9\tau_{0}$')
plt.savefig('jz-multipannel_late_times.pdf',bbox_inches='tight')

plt.show()
'''
'''
ax[1].imshow(jz[1],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[2].imshow(jz[2],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[3].imshow(jz[3],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')

plt.show()

#plt.savefig('jz-mr-2d-map.pdf',bbox_inches='tight')
'''
