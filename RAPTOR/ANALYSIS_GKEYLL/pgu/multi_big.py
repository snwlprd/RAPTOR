import numpy as np
import matplotlib.pyplot as plt
import raptor
import matplotlib.lines as lines



#define empty arrays to fill and key parameters to use 
mr=[25,250,1000,1836]
dd=['init204825','init2048250','init20481000','init20481836']
N_data=250
pi_de = np.zeros((7,2048,2048,1))
pi_di = np.zeros((7,2048,2048,1))
pTh_e = np.zeros((7,2048,2048,1))
pTh_i = np.zeros((7,2048,2048,1))
L = 8*np.pi 
ex, ey = L/2, L/2
alphae = np.zeros(7)
betae = np.zeros(7)
alphai = np.zeros(7)
betai = np.zeros(7)
krts_e = np.zeros(7)
krts_i = np.zeros(7)

for i in range(len(mr)):
    print(i)

    #load in data
    m = mr[i]
    d=dd[i]
    pidel, pThel, pidil, pThil = raptor.pi_d_2d(d,m,204)

    pi_de[i] = pidel
    pi_di[i] = pidil
    pTh_e[i] = pThel
    pTh_i[i] = pThil

    #scale data to make data shown in figure consitant
    alphae[i] = np.max(pi_de[i]) #*0.05)
    alphai[i] = np.max(pi_di[i]) #*0.8)
    
    betae[i] = np.max(pTh_e[i])
    betai[i] = np.max(pTh_i[i])


#font and font size
plt.rcParams.update({'font.size': 9})
#plt.rc('text', usetex=True)


fig,ax=plt.subplots(4,4,sharex='col',sharey='row',figsize=(10,10))
fig.subplots_adjust(wspace=0.025,hspace=0.05)


for j in range(len(mr)):
    
    m=mr[j]
    #pilaat=ax[0][j].imshow(3*pi_de[j]/alpha[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    ax[0][j].imshow(6*pi_de[j]/alphae[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    #ax[0][j].text(-0.9*ex,0.8*ey,str(m))
    ax[0][0].set_ylabel(r'$\Pi^e D^e$')
    ax[0][j].set_xlabel(str(m),size=14) #,bbox=dict(facecolor='none', edgecolor='k', pad=1.0))
    ax[0][j].xaxis.set_label_position('top')
    #ax[j].text(-0.9*ex,-0.8*ey,r'$\kappa=$'+str('%.3g'%krts_e[j]))
    #ax[j].set_xlabel('x $(d_{i})$')

    ax[1][j].imshow(2*pi_di[j]/alphai[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    #ax[1][j].text(-0.9*ex,0.8*ey,str(m))
    ax[1][0].set_ylabel(r'$\Pi^i D^i$')

    ax[2][j].imshow(3*pTh_e[j]/betae[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    #ax[2][j].text(-0.9*ex,0.8*ey,str(m))    
    ax[2][0].set_ylabel(r'$p_e\theta_e$')
    ax[2][j].add_artist(lines.Line2D([0.406*ex,0.406*ex],[-ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[2][j].add_artist(lines.Line2D([-0.102*ex,-0.102*ex],[-ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))
    ax[2][j].add_artist(lines.Line2D([-0.102*ex,0.406*ex],[-ey,-ey],linestyle=':',linewidth=1,color='black'))
    ax[2][j].add_artist(lines.Line2D([-0.102*ex,0.406*ex],[-0.829*ey,-0.829*ey],linewidth=1,linestyle=':',color='black'))

    ra = ax[3][j].imshow(pTh_i[j]/betai[j],cmap='PiYG',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower')
    #ax[3][j].set_xlabel() #('x $(d_{i})$')
    #ax[3][j].xaxis.set_label_position('top') 
    #ax[3][j].text(-0.9*ex,0.8*ey,str(m))
    ax[3][0].set_ylabel(r'$p_i\theta_i$')

    #ax[j][0].set_ylabel('y $(d_{i})$')

    ax[j][0].tick_params(axis='y',direction='in')
    ax[j][0].tick_params(axis='x',direction='in')
    ax[j][0].xaxis.set_ticks_position('both')
    ax[j][0].yaxis.set_ticks_position('both')
    ax[j][0].tick_params(which='minor', direction='in')

    ax[j][1].tick_params(axis='y',direction='in')
    ax[j][1].tick_params(axis='x',direction='in')
    ax[j][1].xaxis.set_ticks_position('both')
    ax[j][1].yaxis.set_ticks_position('both')
    ax[j][1].tick_params(which='minor', direction='in')

    ax[j][2].tick_params(axis='y',direction='in')
    ax[j][2].tick_params(axis='x',direction='in')
    ax[j][2].xaxis.set_ticks_position('both')
    ax[j][2].yaxis.set_ticks_position('both')
    ax[j][2].tick_params(which='minor', direction='in')

    ax[j][3].tick_params(axis='y',direction='in')
    ax[j][3].tick_params(axis='x',direction='in')
    ax[j][3].xaxis.set_ticks_position('both')
    ax[j][3].yaxis.set_ticks_position('both')
    ax[j][3].tick_params(which='minor', direction='in')
#fig.colorbar(pilaat,ax=ax[:],pad=0.02)
#fig.suptitle(r'$p\theta_{electrons}/p\theta_{electrons}^{MAX}$ at $\tau=4.9$')

    #colourbar
fig.colorbar(ra, ax=ax[:],location='right',aspect=50,pad=0.03)#,shrink=1.0)

fig.text(0.44, 0.09, r'x $(d_i)$', ha='center', va='center')

plt.savefig('big-multipannel_v2.pdf',bbox_inches='tight',pad_inches=0.05)

plt.show()
'''
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
'''


