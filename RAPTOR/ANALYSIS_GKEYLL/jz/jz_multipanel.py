import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import raptor
import pandas as pd

mr=[25,250,1000,1836]
N_data=250
lx, ly, lz = 2048, 2048, 1
jz = np.zeros((len(mr),lx,ly,lz))
L = 8*np.pi 
ex, ey = L/2, L/2
krts=np.zeros(len(mr))
dd = ['init204825','init2048250','init20481000','init20481836']
abcd = ['a','b','c','d']

for i in range(len(mr)):
    print(i)

    d = dd[i]
    m = mr[i]
    jxl, jyl, jzl = raptor.j_xyz(d,m,204)

    jz[i] = jzl

    #find kurtosis
    krts[i] = np.mean(jz[i]**4)/(np.mean(jz[i]**2)**2)

#print(np.max(jz))

minmax=[np.min(jz),np.max(jz)]

#pd.read_csv('ENERGIES_2048_' + str(m)+ '_normalised.csv')

fig,ax=plt.subplots(1,4,sharey='row',figsize=(10,3))
fig.subplots_adjust(wspace=0.1)

alpha = np.max(jz)*1.5
print(alpha)

for j in range(len(mr)):
    m=mr[j]
    pilaat=ax[j].imshow(jz[j]/alpha,cmap='seismic',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    ax[j].text(-0.9*ex,0.8*ey,str(m))
    #ax[j].text(0.2,-0.78*ey,'('+str(abcd[j])+')')
    #ax[j].set_xlabel('x $(d_{i})$')
    ax[j].add_artist(lines.Line2D([0.406*ex,0.406*ex],[-ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[j].add_artist(lines.Line2D([-0.102*ex,-0.102*ex],[-ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[j].add_artist(lines.Line2D([-0.102*ex,0.406*ex],[-ey,-ey],linestyle=':',linewidth=1,color='black'))
    ax[j].add_artist(lines.Line2D([-0.102*ex,0.406*ex],[-0.829*ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[j].tick_params(axis='y',direction='in')
    ax[j].tick_params(axis='x',direction='in')
    ax[j].xaxis.set_ticks_position('both')
    ax[j].yaxis.set_ticks_position('both')
    ax[j].tick_params(which='minor', direction='in')


ax[0].set_ylabel('y $(d_{i})$')
# common axis labels
fig.text(0.518, 0.04, r'x $(d_i)$', ha='center', va='center')

fig.colorbar(pilaat,ax=ax[:],location='top',aspect=55,pad=0.04)
#fig.suptitle(r'$j_{z}/(2.28\times10^{-2})$ at $\tau=4.9\tau_{0}$')
plt.savefig('jz-multipannel_late_times.pdf',bbox_inches='tight',pad_inches=0.05)

plt.show()

'''
ax[1].imshow(jz[1],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[2].imshow(jz[2],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[3].imshow(jz[3],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')

plt.show()

#plt.savefig('jz-mr-2d-map.pdf',bbox_inches='tight')
'''
