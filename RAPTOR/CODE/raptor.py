'''
  _____     .       .    _ ___  _______  _____   _____
 |  __ \   /'(  _  )'\  |  __ \|__   __|/     \ |  __ \
 | |__) | / . \/v\/ . \ | |__) |  | |  /   _   \| |__) |
 |  _  / /  _)_`-'_(_  \| .___/   | | |   (_)   |  _  /
 | | \ \/.-~   ).(   ~-.\ |       | |  \       /| | \ \
 |_|  \/'     /\_/\     `\|       |_|   \_____/ |_|  \_\
              "' '"                                    

Rapid Analyisis for Plasma Turbulence Origin Redux    


PYTHON CODE DEFINING FUNCTIONS TO BE USED IN ANALYSIS OF TURBULENT PLASMA SIMULATIONS RAN IN GKEYLL

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
for dd in ['WorkSpace', 'WorkSpace/python', 'AJGAR/TurbAn']:
    if os.path.exists(os.environ['HOME']+'/'+dd):
        sys.path.insert(0,os.environ['HOME']+'/'+dd+'/')
for dd in ['sims/MR_test_closure','sims/MR_test_workload', 'sims/MR_test','sims/MR_sims', 'analysis/MR_analysis', 'analysis/init']:
    if os.path.exists('/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/'):
        sys.path.insert(0,'/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/')
import lpg


def DoubleIntegrate(z):
    '''
        Taking the mean of z
    '''
    total_sum = np.sum(z)
    return total_sum

def pderiv(ar,dx=1,axis=0):
    '''
        Defining the derivative      
    '''
    return (np.roll(ar,1,axis=axis) - np.roll(ar,-1,axis=axis))/(2*dx)

def pdx(ar,dx):
    '''
        Defining dx
    '''
    return pderiv(ar,dx,0)

def pdy(ar,dy):
    '''
        Defining dy
    '''
    return pderiv(ar,dy,1)

def pdz(ar,dz):
    '''
        Defining dz
    '''
    return pderiv(ar,dz,2)

def curl(ax,ay,az,dx,dy,dz):
    '''
        Defining curl
    '''
    cx = pdy(az,dy)-pdz(ay,dz)
    cy = pdz(ax,dz)-pdx(az,dx)
    cz = pdx(ay,dx)-pdy(ax,dy)
    return cx,cy,cz

def curlx(ay,az,dy,dz):
    '''
        Defining curl x
    '''
    return pdy(az,dy)-pdz(ay,dz)

def curly(ax,az,dx,dz):
    '''
        Defining curl y
    '''
    return pdz(ax,dz)-pdx(az,dx)

def curlz(ax,ay,dx,dy):
    '''
        Defining the curl z 
    '''
    return pdx(ay,dx)-pdy(ax,dy)

#def con_alf(x, 

def load_data(mr, j, claws):
    '''
        Loading in the data from the simulations

        mr: specific mass ratio sim
        j: which time step in the sim
        claws: which specific data to load
    '''
    if isinstance(mr, int):
        data = lpg.lpg('initTESTwl'+str(mr))
        data.vars2load(claws)
        data.loadslice(j)
    if isinstance(mr, dict):
        data = {}
        for i in range(len(d)):
            dl = lpg.lpg(str(d[i]))
            dl.vars2load(claws)
            dl.loadslice(j)
            data[i] = {str(d[i]):dl}
    if isinstance(mr, str):
        data = lpg.lpg(mr)
        data.vars2load(claws)
        data.loadslice(j)
    else:
        raise ValueError('Bad value %s (must be int, str or dict)' % mr)
    return data

def load_ddata(d, j, claws):
    '''
        Loading data from simulations with a dictionary, d, and outputting 
        as a dictionary, so that multiple simulations and INIT files can 
        be loaded 
    '''
    data = {}
    
    for i in range(len(d)):
        dl = lpg.lpg(str(d[i]))
        dl.vars2load(claws)
        dl.loadslice(j)
        data[i] = {str(d[i]):dl}
    
    return data
       

def j_xyz(d,mr,j):  #mr, j):
    '''
        Defining the z,y and out of plane, z, current density

        mr: specific mass ratio sim 
        j: which time step in the sim

        NOTE: -mr, mass ratio should be defined elsewhere as mr=mr[m1,mr2,...]
              -CURRENTLY ONLY CALUCLATES OUT OF PLANE CURRENT DENSITY
    ''' 
    #loading in data to be used 
    data = load_data(d,j,('bx','by','bz'))  #load_data(mr, j, ('bx','by','bz'))
    bx, by, bz = data.bx, data.by, data.bz 

    #defining mass ratio, ion inertial length and domain length
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(data.bx)

    #defining grid increment length
    dx = L/gp
    dy = L/gp
    dz = 1e-9

    jx = curlx(by,bz,dy,dz)
    jy = curly(bx,bz,dx,dz)
    jz = curlz(bx,by,dx,dy)

    return jx, jy, jz

def j_dot_E_vxB(d,mr,j):
    '''
    Finding j dot (E + v X B)
    '''
    #loading data 
    data = load_data(d, j, ('ex', 'ey', 'ez','bx','by','bz','jex','jey','jez','rme'))
    ex, ey, ez = data.ex, data.ey, data.ez
    bx, by, bz = data.bx, data.by, data.bz
    jex, jey, jez =data.jex, data.jey, data.jez
    rme = data.rme
    
    jx, jy, jz = j_xyz(mr,j)

    #defining mass ratio, ion inertial length and domain length
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(data.bx)

    #defining grid increment length
    dx = L/gp
    dy = L/gp
    dz = 1e-9

    #finding v
    vx, vy, vz = jex/rme, jey/rme, jez/rme

    #v X B
    #jvBx = jx*((vy*bz)-(vz*by))
    #jvBy = -1*vx*((jy*bz)-(jz*by))
    #jvBz = bx*((jx*vz)-(jz*vy))
    
    jvB = jx*((vy*bz)-(vz*by)) - vx*((jy*bz)-(jz*by)) + bx*((jx*vz)-(jz*vy))
    #jvB = jx*((vy*bz)-(vz*by)) - jy*((vx*bz)-(vz*bx)) + jz*((vx*by)-(vy*bx))

    #j dot E
    j_dot_E = jx*ex + jy*ey + jz*ez

    #j dot (E + v X B)
    jEvB = j_dot_E + jvB

    return jEvB

def j_dot_E(d,mr, j):
    '''
    Finding j dot E
    '''
    data = load_data(d, j, ('ex', 'ey', 'ez','bx','by','bz'))
    ex, ey, ez = data.ex, data.ey, data.ez
    bx, by, bz = data.bx, data.by, data.bz

    #defining mass ratio, ion inertial length and domain length
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(bx)

    #defining grid increment length
    dx = L/gp
    dy = L/gp
    dz = 1e-9

    jx, jy, jz = j_xyz(d,mr,j)

    j_dot_E = jx*ex + jy*ey + jz*ez

    return j_dot_E

def vort_xyz(d,mr,j):
    '''
        Defining a function to work out vorticity
    '''
    #loading in data to be used 
    data = load_data(d, j, ('jex','jey','jez','jix','jiy','jiz','rme','rmi'))
    jex, jey, jez = data.jex, data.jey, data.jez
    jix, jiy, jiz = data.jix, data.jiy, data.jiz
    rme, rmi = data.rme, data.rmi

    #defining mass ratio, ion inertial length and domain length
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(jex)

    #defining grid increment
    dx = L/gp
    dy = L/gp

    #defining components of the flow
    uex = jex/rme
    uey = jey/rme
    uez = jez/rme

    uix = jix/rmi
    uiy = jiy/rmi
    uiz = jiz/rmi

    #calculating electron (vez) and ion (viz) vorticity 
    vez = curlz(uex,uey,dx,dy)
    viz = curlz(uix,uiy,dx,dy)

    return vez, viz

def tci_tnl(d,mr,j): #,l):
    '''
        Defining a function to calculate the ratio of ion cyclotron timescale 
        to non linear timescale, at a length l

        mr: specific mass ratio sim 
        j: which time step in the sim
        l: length scale

        NOTE: -mr, mass ratio should be defined elsewhere as mr=mr[m1,mr2,...]
              -ion-cyclotron timescale, tci, is the inverse of the 
               ion-cylotron frequency
              -ONLY WORKING OUT SCALE d_i CURRENTLY
    '''
    #ion intertial length
    d_i=np.sqrt(mr)
    
    #loading in data to be used 
    data = load_data(d, j, ('bx','by','bz'))
    bx, by, bz = data.bx, data.by, data.bz 
    #bxs, bys, bzs = (bx/np.sqrt(4*np.pi*mr)), (by/np.sqrt(4*np.pi*mr)), (bz/np.sqrt(4*np.pi*mr))
  
    #calculating magnetic field magnitude
    Bmag = np.sqrt(bx**2 + by**2 + bz**2)

    #calculating ion-cylotron timescale
    tci = mr/(Bmag)
    '''
        #calculating magnetic spectrum
        k,ebx,eby,ebz,eb=af.PerpSpecVec(bx,by,bz,lx,ly,lz)
        # Equivalent wavenumber
        k_l = 1./l
        # Energy at that wavenumber
        idxk_l= np.argmin(np.abs(k-k_l))
        ek_l=eb[idxk_l]
        #calculating non-linear timescale at l
        tnl = 1/(k_l*np.sqrt(ek_l))

        return tc/tnl, tci, tnl
    '''
    #calculating dB
    dBx = np.roll(bx,22,axis=0)-bx
    dBy = np.roll(by,22,axis=0)-by
    dBz = np.roll(bz,22,axis=0)-bz

    dB = np.sqrt(dBx**2+dBy**2+dBz**2)
    dv_A = dB/(np.sqrt(4*np.pi*mr))

    #calculating non-linear timescale at d_i
    tnl = d_i/dv_A

    #calculating ratio of ion-cyclotron to non-linear
    tci_tnl = tci/tnl
                                        
    return tci, tnl, tci_tnl


def qpq(mr):
    '''
        Defining a function to calculate the total heating, qi+qe
    '''
    #loading in data from chosen file
    f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/comp_other/energies/ENERGIES_2048_' + str(mr) + '.csv')
    
    # Average derivative value for all the energies
    a=f.diff(periods=1,axis=0).mean()

    #calculating qi+qe
    qpq = a['Eith'] + a['Eeth']
    
    return qpq

def qi_qe(mr): #,j):
    '''
        Defining a function to calculate the time rate of change 
        of ion to electron internal energy, qi/qe
        
        mr: specific mass ratio sim 
        j: which time step in the sim

        'NOTE: -mr, mass ratio should be defined elsewhere, mr=[m1,mr2,...]'
    '''    
    #loading in data from chosen file
    f = pd.read_csv('/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_test/analysis/comp_other/energies/ENERGIES_2048_' + str(mr) + '.csv')
   
    #calculating tau_naught
    tau_zero = (8*np.pi*np.sqrt(mr))/(2*np.pi*0.0023)

    #calculating dt
    dt = 6*tau_zero/250

    #Average derivative value for all the energies
    #a=f[167:250].diff(periods=1,axis=0).mean()
    #a=f.rolling(10).mean().diff(periods=1,axis=0)/dt
    a=f.diff(periods=1,axis=0)/dt #(where dt is the time between each time slice
    qi = a['Eithn']
    qe = a['Eethn']
    qi_qe = (a['Eithn'])/(a['Eethn']) 
    qpq = a['Eithn'] + a['Eethn']
    
    #N_data = 250
    #eEeth =[]
    #eEith = []
    #for j in range(N_data - 1):
    #    Eeth = f['Eeth'][j+1] - f['Eeth'][j]
    #    Eith = f['Eith'][j+1] - f['Eith'][j]    
    #    eEeth.append(Eeth)
    #    eEith.append(Eith)
    #aEeth = np.mean(eEeth)
    #aEith = np.mean(eEith) 
    #qi_qe = aEith/aEeth
    #qpq = aEith + aEeth

    return qe, qi, qi_qe, qpq  

def p_g_u_2d(d,mr,j):
    '''
        Defining a function to calculate P dot Grad dot u with 2 
        spatial dimentions, x,y
       
        mr: specific mass ratio sim 
        j: which time step in the sim

        NOTE: -mr, mass ratio should be defined elsewhere as mr=mr[m1,mr2,...]
              -gp, number of grid points should be defined elsewhere as gp=...
    ''' 
    #loading in data to be used 
    data = load_data(d, j, ('jex', 'jey', 'jez', 'jix', 'jiy', 'jiz', 'pexx', 'peyy', 'pezz', 'pixx', 'piyy', 'pizz','pixy','pixz','piyz', 'pexy','pexz','peyz','rme', 'rmi'))

    #defining mass ratio, ion inertial length and domain length        
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(data.jex)

    #defining grid step size in x,y and z direction
    dx=L/gp
    dy=L/gp
    dz=1e-6

    #defining each component of the symetric pressure tensor
    pexx, pexy, pexz = data.pexx, data.pexy, data.pexz 
    peyx, peyy, peyz = data.pexy, data.peyy, data.peyz
    pezx, pezy, pezz = data.pexz, data.peyz, data.pezz  

    pixx, pixy, pixz = data.pixx, data.pixy, data.pixz
    piyx, piyy, piyz = data.pixy, data.piyy, data.piyz
    pizx, pizy, pizz = data.pixz, data.piyz, data.pizz

    #calculating u
    uex = (data.jex)/data.rme
    uey = (data.jey)/data.rme
    uez = (data.jez)/data.rme

    uix = (data.jix)/data.rmi
    uiy = (data.jiy)/data.rmi
    uiz = (data.jiz)/data.rmi

    #calculating P dot Grad dot u for ions and electrons
    epgu = -1*(  (pexx*pdx(uex,dx)) + (pexy*pdx(uey,dx)) + (pexz*pdx(uez,dx)) \
               + (peyx*pdy(uex,dy)) + (peyy*pdy(uey,dy)) + (peyz*pdy(uez,dy))) 
                
    ipgu = -1*(  (pixx*pdx(uix,dx)) + (pixy*pdx(uiy,dx)) + (pixz*pdx(uiz,dx)) \
               + (piyx*pdy(uix,dy)) + (piyy*pdy(uiy,dy)) + (piyz*pdy(uiz,dy))) 
                                                     
    #calculating the combined ion and electron P dot Grad dot u 
    apgu = epgu+ipgu
    
    return epgu, ipgu, apgu

def p_g_u_3d(d,mr,j):
    '''
        Defining a function to calculate P dot Grad dot u with 3 spatial dimentions
            mr: specific mass ratio sim                
            j: is at which time step you are at in the sim

            NOTE: -mr, mass ratio should be defined else where as mr=mr[m1,mr2,...]
                  -gp, number of grid points should be defined elsewhere as gp=...
    '''
    #loading in data to be used 
    data = load_data(d, j, ('jex', 'jey', 'jez', 'jix', 'jiy', 'jiz', 'pexx', 'peyy', 'pezz', 'pixx', 'piyy', 'pizz','pixy','pixz','piyz', 'pexy','pexz','peyz','rme', 'rmi'))

    #defining mass ratio, ion inertial length and domain length        
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(data.jex)

    #defining grid step size in x,y and z direction
    dx=L/gp
    dy=L/gp
    dz=1e-6

    #defining each component of the symetric pressure tensor
    pexx, pexy, pexz = data.pexx, data.pexy, data.pexz
    peyx, peyy, peyz = data.pexy, data.peyy, data.peyz
    pezx, pezy, pezz = data.pexz, data.peyz, data.pezz

    pixx, pixy, pixz = data.pixx, data.pixy, data.pixz
    piyx, piyy, piyz = data.pixy, data.piyy, data.piyz
    pizx, pizy, pizz = data.pixz, data.piyz, data.pizz

    #calculating u
    uex = (data.jex)/data.rme
    uey = (data.jey)/data.rme
    uez = (data.jez)/data.rme
    
    uix = (data.jix)/data.rmi
    uiy = (data.jiy)/data.rmi
    uiz = (data.jiz)/data.rmi

    #calculating P dot Grad dot u for ions and electrons
    epgu = -1*(  (pexx*pdx(uex,dx)) + (pexy*pdx(uey,dx)) + (pexz*pdx(uez,dx)) \
               + (peyx*pdy(uex,dy)) + (peyy*pdy(uey,dy)) + (peyz*pdy(uez,dy)) \
               + (pezx*pdz(uex,dz)) + (pezy*pdz(uey,dy)) + (pezz*pdz(uez,dz))  )


    ipgu = -1*(  (pixx*pdx(uix,dx)) + (pixy*pdx(uiy,dx)) + (pixz*pdx(uiz,dx)) \
               + (piyx*pdy(uix,dy)) + (piyy*pdy(uiy,dy)) + (piyz*pdy(uiz,dy)) \
               + (pizx*pdz(uix,dz)) + (pizy*pdz(uiy,dy)) + (pizz*pdz(uiz,dz))  )
    
    #calculating the combined ion and electron P dot Grad dot u 
    apgu = epgu+ipgu

    return epgu, ipgu, apgu

def pi_d_2d(d,mr,j):
    '''
        Defining a function to calculate Pi-D with 2 spatial
        dimentions, x,y
                               
        mr: specific mass ratio sim 
        j: which time step in the sim

        NOTE: -mr, mass ratio should be defined elsewhere as mr=mr[m1,mr2,...]
    '''
    #loading in data to be used 
    data = load_data(d, j, ('jex', 'jey', 'jez', 'jix', 'jiy', 'jiz', 'pexx', 'peyy', 'pezz', 'pixx', 'piyy', 'pizz','pixy','pixz','piyz', 'pexy','pexz','peyz','rme', 'rmi'))

    #defining mass ratio, ion inertial length and domain length        
    d_i = np.sqrt(mr)
    L = 8.*np.pi*d_i
    gp = len(data.jex)

    #defining grid step size in x,y and z direction
    dx=L/gp
    dy=L/gp
    dz=1e-6

    #defining each component of the symetric pressure tensor
    pexx, pexy, pexz = data.pexx, data.pexy, data.pexz
    peyx, peyy, peyz = data.pexy, data.peyy, data.peyz
    pezx, pezy, pezz = data.pexz, data.peyz, data.pezz
    
    pixx, pixy, pixz = data.pixx, data.pixy, data.pixz
    piyx, piyy, piyz = data.pixy, data.piyy, data.piyz
    pizx, pizy, pizz = data.pixz, data.piyz, data.pizz

    #calculating u
    uex = (data.jex)/data.rme
    uey = (data.jey)/data.rme
    uez = (data.jez)/data.rme

    uix = (data.jix)/data.rmi
    uiy = (data.jiy)/data.rmi
    uiz = (data.jiz)/data.rmi

    
    #calculating scalar pressure
    spe = (1/3)*(pexx+peyy+pezz)
    spi = (1/3)*(pixx+piyy+pizz)

    #calculating the deviatoric pressure tensor
    Piexx, Piexy, Piexz = pexx - spe, pexy, pexz
    Pieyx, Pieyy, Pieyz = peyx, peyy - spe, peyz
    Piezx, Piezy, Piezz = pezx, pezy, pezz - spe

    Piixx, Piixy, Piixz = pixx - spi, pixy, pixz
    Piiyx, Piiyy, Piiyz = piyx, piyy - spi, piyz
    Piizx, Piizy, Piizz = pizx, pizy, pizz - spi

    #calculating strain-rate tensor
    #Sexx, Sexy, Sexz = (1/2)*(pdx(uex,dx)+pdx(uex,dx)), (1/2)*(pdx(uey,dx)+pdy(uex,dy)), (1/2)*(pdx(uez,dx)+pdz(uex,dz))
    #Seyx, Seyy, Seyz = (1/2)*(pdy(uex,dy)+pdx(uey,dx)), (1/2)*(pdy(uey,dy)+pdy(uey,dy)), (1/2)*(pdy(uez,dy)+pdz(uey,dz))
    #Sezx, Sezy, Sezz = (1/2)*(pdz(uex,dz)+pdx(uez,dx)), (1/2)*(pdz(uey,dz)+pdy(uez,dy)), (1/2)*(pdz(uez,dz)+pdz(uez,dz))

    Sexx, Sexy, Sexz = pdx(uex,dx), (1/2)*(pdx(uey,dx)+pdy(uex,dy)), (1/2)*(pdx(uez,dx))
    Seyx, Seyy, Seyz = (1/2)*(pdy(uex,dy)+pdx(uey,dx)), pdy(uey,dy), (1/2)*(pdy(uez,dy))
    Sezx, Sezy, Sezz = (1/2)*(pdx(uez,dx)), (1/2)*(pdy(uez,dy)), 0

    Sixx, Sixy, Sixz = pdx(uix,dx), (1/2)*(pdx(uiy,dx)+pdy(uix,dy)), (1/2)*(pdx(uiz,dx))
    Siyx, Siyy, Siyz = (1/2)*(pdy(uix,dy)+pdx(uiy,dx)), pdy(uiy,dy), (1/2)*(pdy(uiz,dy))
    Sizx, Sizy, Sizz = (1/2)*(pdx(uiz,dx)), (1/2)*(pdy(uiz,dy)), 0

    #calculating dilation, Theta
    The = Sexx + Seyy + Sezz #Th = pdx(uex,dx)+pdy(uey,dy) #+pdz(uez,dz)
    Thi = Sixx + Siyy + Sizz

    #calculating the traceless strain-rate tensor
    Dexx, Dexy, Dexz = Sexx - (1/3)*The, Sexy, Sexz
    Deyx, Deyy, Deyz = Seyx, Seyy - (1/3)*The, Seyz
    Dezx, Dezy, Dezz = Sezx, Sezy, Sezz - (1/3)*The

    Dixx, Dixy, Dixz = Sixx - (1/3)*Thi, Sixy, Sixz
    Diyx, Diyy, Diyz = Siyx, Siyy - (1/3)*Thi, Siyz
    Dizx, Dizy, Dizz = Sizx, Sizy, Sizz - (1/3)*Thi
            

    #calculating p Theta
    pTh_e = -1*(spe*The)
    pTh_i = -1*(spi*Thi)

    #calculating Pi-D
    pi_d_e = -1*((Piexx*Dexx) + (Piexy*Dexy) + (Piexz*Dexz) + (Pieyx*Deyx) + (Pieyy*Deyy) + (Pieyz*Deyz) + (Piezx*Dezx) + (Piezy*Dezy) + (Piezz*Dezz)) 
    pi_d_i = -1*((Piixx*Dixx) + (Piixy*Dixy) + (Piixz*Dixz) + (Piiyx*Diyx) + (Piiyy*Diyy) + (Piiyz*Diyz) + (Piizx*Dizx) + (Piizy*Dizy) + (Piizz*Dizz))

    return pi_d_e, pTh_e, pi_d_i, pTh_i

















