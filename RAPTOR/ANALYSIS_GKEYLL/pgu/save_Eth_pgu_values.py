import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import raptor

mr = [25,50,100,250,500,1000,1836]
dd = ['init204825','init204850','init2048100','init2048250','init2048500','init20481000','init20481836']

N_data = 250 

t = np.linspace(0,6,250)

for i in range(len(mr)):
    print(i)
    d = dd[i]
    m = mr[i]
    
    #setting up list to append j dot E(jde), sun of j dot E(sjde), PiD(pid), pTheta(pth), sum of PiD(spid) and sum of Pth(spth) foe ions(i) and elctrons(e)
    pidi, pth_i,  pide, pth_e = [], [], [], []
    spidi, spide, spth_i, spth_e = [], [], [], []
    jde, sjde = [], []

    #setting up p dot grad dot u as above
    pgui, pgue = [], []
    spgui, spgue = [], []

    #loading in Energy data for mr sim
    f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_'+str(m)+'_normalised.csv')
    g = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_'+str(m)+'.csv')

    for j in range(0,N_data):
        print(j)
        
        #getting PiD(pi_d, pd, pid), and pTheta(pth) for ions(i) and elecrtons(e)
        pi_de, pthe, pi_di, pthi = raptor.pi_d_2d(d,m,j)
        pde_a, pthe_a, pdi_a, pthi_a = np.mean(pi_de), np.mean(pthe), np.mean(pi_di), np.mean(pthi)
        pidi.append(pdi_a), pth_i.append(pthi_a), pide.append(pde_a), pth_e.append(pthe_a)

        #getting P dot grad dot u for mr sim
        et, it, at = raptor.p_g_u_2d(d,m,j)
        ita, eta = np.mean(it), np.mean(et)
        pgui.append(ita), pgue.append(eta)

        #getting j dot E
        l_jde = raptor.j_dot_E(d,m,j)
        a_jde = np.mean(abs(l_jde))*(((4*np.pi/m)**(3/2))*m)/np.sqrt(m)
        jde.append(a_jde)

        #creating list of cumulative average PiD and pTheta
        if j ==0:
            spidi.append(pidi[j]), spide.append(pide[j]), spth_i.append(pth_i[j]), spth_e.append(pth_e[j])
            spgui.append(pgui[j]), spgue.append(pgue[j])        
            sjde.append(jde[j])
        else:
            sum_pidi = pidi[j] + spidi[j-1]
            sum_pide = pide[j] + spide[j-1]
            sum_pthi = pth_i[j] + spth_i[j-1]
            sum_pthe = pth_e[j] + spth_e[j-1]
            sum_pgui = pgui[j] + spgui[j-1]
            sum_pgue = pgue[j] + spgue[j-1]
            sum_jde = jde[j] + sjde[j-1]

            spidi.append(sum_pidi), spide.append(sum_pide), spth_i.append(sum_pthi), spth_e.append(sum_pthe)
            spgui.append(sum_pgui), spgue.append(sum_pgue)
            sjde.append(sum_jde)

    #calculating tau_naught
    tau_zero = (8*np.pi*np.sqrt(m))/(2*np.pi*0.0023)
    
    #calculating dt
    dt = 6*tau_zero/250 
    
    #turning list into an array
    spidi_np, spide_np, spthi_np, spthe_np = np.array(spidi), np.array(spide), np.array(spth_i), np.array(spth_e)
    spgui_np, spgue_np = np.array(spgui), np.array(spgue)
    sjde_np = np.array(sjde)

    #cumulative P dot grad dot u with dt
    c_PiD_ions = -1*dt*spidi_np
    c_PiD_elec = -1*dt*spide_np
    c_pTh_ions = -1*dt*spthi_np
    c_pTh_elec = -1*dt*spthe_np

    #cumulative P dot grad dot u with dt
    c_pgu_ions = -1*dt*spgui_np
    c_pgu_elec = -1*dt*spgue_np

    #value to normlaise
    norm = g['Eb'][0] + g['Eefl'][0] + g['Eifl'][0]

    #save values
    V  = {'PiD_ions':c_PiD_ions, 'PiD_elec':c_PiD_elec, 'pTh_ions':c_pTh_ions, 'pTh_elec':c_pTh_elec, 'pgu_ions':c_pgu_ions, 'pgu_elec':c_pgu_elec, 'norm':norm}
    df = pd.DataFrame(V)

    df.to_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/pgu/PI_D_VALUES_{}.csv'.format(str(m)))


'''
    plt.plot(t,(c_PiD_ions)/(norm),label=r'$\Pi D$')
    plt.plot(t,(c_pTh_ions)/(norm),label=r'$p \Theta$')
    plt.plot(t,(c_pgu_ions)/(norm),label=r'$(P \cdot \nabla) \cdot u$')
    #plt.plot(t,sjde_np/(f['Eith'][249]-f['Eith'][0]),label=r'$j \cdot E$')
    plt.plot(t,(g['Eith']-g['Eith'][0])/(norm),linestyle='--',label=r'$\Delta \epsilon^{th}_{i}$')
    plt.legend()
    plt.title(r'Ions with $m_{i}/m_{e}$ = ' + str(m))
    plt.xlabel(r'Eddy turnover time, $\tau_{0}$')
    plt.ylabel(r'Cumulative $\Pi D$, Cumulative $p \Theta$')
    plt.savefig('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_test/analysis/pid/cumulative_PiD_pTh_vs_Eith'+str(m)+'.pdf')
#    plt.show()
    plt.clf()
    
    plt.plot(t,(c_PiD_elec)/(norm),label=r'$\Pi D$')
    plt.plot(t,(c_pTh_elec)/(norm),label=r'$p \Theta$')
    plt.plot(t,(c_pgu_elec)/(norm),label=r'$(P \cdot \nabla) \cdot u$')
    plt.plot(t,(sjde_np)/(norm),linestyle='-.',label=r'$j \cdot E$')
    plt.plot(t,(g['Eeth']-g['Eeth'][0])/(norm),linestyle='--',label=r'$\Delta \epsilon^{th}_{e}$')
    plt.legend()
    plt.title(r'Electrons with $m_{i}/m_{e}$ = ' + str(m))
    plt.xlabel(r'Eddy turnover time, $\tau_{0}$')
    plt.ylabel(r'Cumulative $\Pi D$, Cumulative $p \Theta$')
    plt.savefig('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_test/analysis/pid/cumulative_PiD_pTh_vs_Eeth'+str(m)+'.pdf')
#    plt.show()
    plt.clf()
'''
