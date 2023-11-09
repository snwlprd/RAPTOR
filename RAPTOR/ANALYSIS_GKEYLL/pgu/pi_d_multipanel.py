import numpy as np
import matplotlib.pyplot as plt
import raptor

#define empty arrays to fill and key parameters to use 
mr=[25,250,1000,1836]
dd=['init204825','init2048250','init20481000','init20481836']
N_data=250
pi_de = np.zeros((7,2048,2048,1))
pi_di = np.zeros((7,2048,2048,1))
L = 8*np.pi 
ex, ey = L/2, L/2
alpha = np.zeros(7)
beta = np.zeros(7)
krts_e = np.zeros(7)
krts_i = np.zeros(7)

for i in range(len(mr)):
    print(i)

    #load in data
    m = mr[i]
    d=dd[i]
    pidel, pThel, pidil, pThil = raptor.pi_d_2d(d,m,204)

    pi_de[i] = pThel  #pidel
    pi_di[i] = pThil  #pidil

    #scale data to make data shown in figure consitant
    alpha[i] = np.max(pi_de[i]) #*0.05)
    beta[i] = np.max(pi_di[i]) #*0.8)

    #find kurtosis
    krts_e[i] = np.mean((pi_de[i]/alpha[i])**4)/(np.mean((pi_de[i]/alpha[i])**2)**2)
    krts_i[i] = np.mean((pi_di[i]/beta[i])**4)/(np.mean((pi_di[i]/beta[i])**2)**2)        

print(krts_e)
print(krts_i)

fig,ax=plt.subplots(1,4,sharey='row',figsize=(10,2))

for j in range(len(mr)):
    m=mr[j]
    pilaat=ax[j].imshow(3*pi_de[j]/alpha[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    ax[j].text(-0.9*ex,0.8*ey,str(m))
    ax[j].text(-0.9*ex,-0.8*ey,r'$\kappa=$'+str('%.3g'%krts_e[j]))
    ax[j].set_xlabel('x $(d_{i})$')



ax[0].set_ylabel('y $(d_{i})$')
fig.colorbar(pilaat,ax=ax[:],pad=0.02)
fig.suptitle(r'$p\theta_{electrons}/p\theta_{electrons}^{MAX}$ at $\tau=4.9$')
plt.savefig('pTh_e-multipannel.pdf',bbox_inches='tight')

plt.show()

plt.clf()

fig,ax=plt.subplots(1,4,sharey='row',figsize=(10,2))

for j in range(len(mr)):
    m=mr[j]
    pilaat=ax[j].imshow(2*pi_di[j]/beta[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    ax[j].text(-0.9*ex,0.8*ey,str(m))
    ax[j].text(-0.9*ex,-0.8*ey,r'$\kappa=$'+str('%.3g'%krts_i[j]))
    ax[j].set_xlabel('x $(d_{i})$')

ax[0].set_ylabel('y $(d_{i})$')
fig.colorbar(pilaat,ax=ax[:],pad=0.02)
fig.suptitle(r'$p\theta_{ions}/p\theta_{ions}^{MAX} $ at $\tau=4.9$')  #5.24$')
plt.savefig('pTh_i-multipannel.pdf',bbox_inches='tight')

plt.show()



