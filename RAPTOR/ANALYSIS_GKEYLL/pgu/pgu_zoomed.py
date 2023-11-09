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
pth = np.zeros((len(mr),lx,ly,lz))
L = 8*np.pi 
ex, ey = L/2, L/2
krts=np.zeros(len(mr))
dd = ['init204825','init2048250','init20481000','init20481836'] 
abcd = ['a','b','c','d']
alpha = np.zeros(7)

for i in range(len(mr)):
    print(i)

    d = dd[i]
    m = mr[i]
    pidel, pThel, pidil, pThil = raptor.pi_d_2d(d,m,204)

    pth[i] = pThel

    alpha[i] = np.max(pth[i])
    #find kurtosis
    #krts[i] = np.mean(jz[i]**4)/(np.mean(jz[i]**2)**2)

#print(np.max(jz))

minmax=[np.min(pth),np.max(pth)]

#font and font size
plt.rcParams.update({'font.size': 9})


fig, ax = plt.subplots(4,1,sharex='col',figsize=(3.5,5))
fig.subplots_adjust(hspace=0.05)


#alpha = np.max(pth)
print(alpha)

#plt.imshow(jz[0,0:175,920:1440]/alpha,cmap='seismic',vmin=-1,vmax=1,origin='lower')
#plt.savefig('jz-zoomed_1000_late_times.pdf',bbox_inches='tight')

#plt.show()

#fig,ax=plt.subplots(1,4,sharey='row',figsize=(10,2))

for j in range(len(mr)):
    m=mr[j]  #920:1440
    
    #pilaat=ax[j].imshow(1.5*pth[j][0:175,920:1440]/alpha[j],cmap='PiYG',interpolation='none',vmin=-1,vmax=1,origin='lower') #,extent=[-ex,ex,-ey,ey],origin='lower')
    pilaat=ax[j].imshow(3*pth[j][0:175,920:1440]/alpha[j],cmap='PiYG',interpolation='none',vmin=-1,vmax=1,origin='lower')
    ax[j].text(15,140,str(m))  #512:1536  #'('+str(abcd[j])+')') #str(m))
    
    #ax[j].text(410,110,r'('+str(abcd[j])+')') #+str('%.3g'%krts[j]))
#    ax[j].set_xlabel('x $(d_{i})$')

#ax[0].set_ylabel('y $(d_{i})$')
#fig.colorbar(pilaat,ax=ax[:],pad=0.02)
#fig.suptitle(r'$j_{z}/(2.45\times10^{-2})$ at $\tau=4.9\tau_{0}$')
for axx in ax:
    axx.set_yticks([])
    axx.set_xticks([])

plt.savefig('pth-zoomed_middle_multipannel_late_times.pdf',bbox_inches='tight') #, pad_inches=0)

plt.show()

'''
ax[1].imshow(jz[1],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[2].imshow(jz[2],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')
ax[3].imshow(jz[3],cmap='seismic',vmin=minmax[0],vmax=minmax[1],extent=[-ex,ex,-ey,ey],origin='lower')

plt.show()

#plt.savefig('jz-mr-2d-map.pdf',bbox_inches='tight')
'''
