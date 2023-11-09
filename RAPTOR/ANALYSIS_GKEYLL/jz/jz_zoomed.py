import numpy as np
#import matplotlib as mpl
#mpl.rc('text', usetex=True)
#mpl.rc('font', family='serif', serif='cm10', size=10)
#import matplotlib.lines as lines
import matplotlib.pyplot as plt
import raptor

mr=[25,250,1000,1836] #, 1836] 
N_data=250
lx, ly, lz = 2048, 2048, 1 #[550,2048], [550,2048], 1
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

#font and font size
plt.rcParams.update({'font.size': 9})
#plt.rc('text', usetex=True)

fig, ax = plt.subplots(4,1,sharex='col',figsize=(3.5,5))
fig.subplots_adjust(hspace=0.05)


alpha = np.max(jz)*1.5
print(alpha)

#plt.imshow(jz[0,0:175,920:1440]/alpha,cmap='seismic',vmin=-1,vmax=1,origin='lower')
#plt.savefig('jz-zoomed_1000_late_times.pdf',bbox_inches='tight')

#plt.show()

#fig,ax=plt.subplots(1,4,sharey='row',figsize=(10,2))

for j in range(len(mr)):
    m=mr[j]
    pilaat=ax[j].imshow(jz[j][0:175,920:1440]/alpha,cmap='seismic',interpolation='none',vmin=-1,vmax=1,origin='lower') #,extent=[-ex,ex,-ey,ey],origin='lower')
    ax[j].text(15,140,str(m))   #'('+str(abcd[j])+')') #str(m))
    #ax[j].text(410,110,r'('+str(abcd[j])+')') #+str('%.3g'%krts[j]))
#    ax[j].set_xlabel('x $(d_{i})$')

for axx in ax:
    axx.set_yticks([])
    axx.set_xticks([])


#ax[0].set_ylabel('y $(d_{i})$')
#fig.colorbar(pilaat,ax=ax[:],pad=0.02)
#fig.suptitle(r'$j_{z}/(2.45\times10^{-2})$ at $\tau=4.9\tau_{0}$')

#set ticks inward
#ax[0].tick_params(axis='y',direction='in')
#ax[0].tick_params(axis='x',direction='in')
#ax[0].xaxis.set_ticks_position('both')
#ax[0].yaxis.set_ticks_position('both')
#ax[0].tick_params(which='minor', direction='in')

#ax[1].tick_params(axis='y',direction='in')
#ax[1].tick_params(axis='x',direction='in')
#ax[1].xaxis.set_ticks_position('both')
#ax[1].yaxis.set_ticks_position('both')
#ax[1].tick_params(which='minor', direction='in')

#ax[2].tick_params(axis='y',direction='in')
#ax[2].tick_params(axis='x',direction='in')
#ax[2].xaxis.set_ticks_position('both')
#ax[2].yaxis.set_ticks_position('both')
#ax[2].tick_params(which='minor', direction='in')

#ax[3].tick_params(axis='y',direction='in')
#ax[3].tick_params(axis='x',direction='in')
#ax[3].xaxis.set_ticks_position('both')
#ax[3].yaxis.set_ticks_position('both')
#ax[3].tick_params(which='minor', direction='in')


plt.savefig('jz-zoomed_multipannel_late_times.pdf',bbox_inches='tight') #, pad_inches=0)

plt.show()

'''
ax[1].imshow(jz[1],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[2].imshow(jz[2],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[3].imshow(jz[3],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')

plt.show()

#plt.savefig('jz-mr-2d-map.pdf',bbox_inches='tight')
'''
