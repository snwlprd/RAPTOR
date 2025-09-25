'''
               .       .                                  
              /'(  _  )'\                       
             / . \/v\/ . \                     
            /  _)_`-'_(_  \                        
           /.-~   ).(   ~-.\                   
          /'     /\_/\     `\                   
              /| "' '" |\                       
             /||| ||| |||\                     
            /\/\/\/\/\/\/<^>                   
           /        <^>   |\                   
         <^>         |      \                    
         /|   <^>            \                 
       <^>     |              <^>              
    <^><^><^><^><^><^><^><^><^><^><^>            
  <^><^><^><^><^><^><^><^><^><^><^><^>                       
 |  ____| \   / |  __ \|__   __|  ____|
 | |___  \ \ / /| |__) |  | |  | |___
 |  ___|  \   / |  _  /   | |  |  ___|
 | |____   | |  |   \ \ __| |__| |____
 |______|  |_|  |_|  \_\_______|______|           
 
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage

#derivatives
def pderiv(ar,dx=1,axis=0):
    '''
        Defining the derivative      
    '''
    if dx==0:
        pderiv = 0
    else:
        pderiv = (np.roll(ar,1,axis=axis) - np.roll(ar,-1,axis=axis))/(2*dx)
    return pderiv

def pdx(ar,dx):
    '''
        Calculating dx
    '''
    return pderiv(ar,dx,0)

def pdy(ar,dy):
    '''
        Calculating dy
    '''
    return pderiv(ar,dy,1)

def pdz(ar,dz):
    '''
        Calculating dz
    '''
    return pderiv(ar,dz,2)

#curls
def curl(ax,ay,az,dx,dy,dz):
    '''
        Calculating curl
    '''
    cx, cy, cz = pdy(az,dy)-pdz(ay,dz), pdz(ax,dz)-pdx(az,dx), pdx(ay,dx)-pdy(ax,dy)
    
    return cx,cy,cz
    
#kurtosis
def kts(d):
    '''
        Calculating kurtosis
    '''
    krt = np.zeros(len(d)) 
    for i in range(len(d)):
        k = np.mean(d[i]**4)/np.mean(d[i]**2)**2
        krt[i] = k

    return krt

#Germano and Favre filtering
def germano(f,l):
    '''
    Fucntion to perform Germano Filtering with a gaussian
    l-filtered scale
    ar-quantity to filter
    '''
    gf_f = ndimage.gaussian_filter(f[:,:,0],l,mode='wrap')

    return gf_f

def favre(f,rho,l):
    '''
    Function to perform Favre Filtering
    '''
    a, b = germano(f*rho,l), germano(rho,l)
    
    return a/b
    
#define velocity
def u(px,py,pz,rm):
    '''
        Calculating velocity, v = (m*n*v)/(m*n) 
    '''
    ux, uy, uz = px/rm, py/rm, pz/rm
    
    return ux, uy, uz

#define vorticity, v = Nable X u 
def vort(px,py,pz,rm,dx,dy,dz=0):
    '''
        Calculating Vorticity
    '''
    ux, uy, uz = u(px,py,pz,rm)
    vx, vy, vz = curlx(uy,uz,dy,dz), curly(ux,uz,dx,dz), curlz(ux,uy,dx,dy)

    return vx, vy, vz

#current
def j(bx,by,bz,dx,dy,dz=0):
    '''
        Calculating Current Density, j = nabla X B
    '''
    jx, jy, jz = curlx(by,bz,dy,dz), curly(bx,bz,dx,dz), curlz(bx,by,dx,dy)

    return jx, jy, jz
    
def jdE(ex,ey,ez,bx,by,bz,dx,dy,dz=0):
    '''
        Calculating j dot E
    '''	
    jx, jy, jz = j(bx,by,bz,dx,dy,dz)
    j_d_E = jx*ex + jy*ey + jz*ez

    return j_d_E

#energies 
def EB(bx,by,bz,eps=1): 
    '''
        Calculating Magnetic Field Energy 
    '''
    E_B = 0.5*eps*np.mean(bx**2 + by**2 + bz**2)
    
    return E_B

def EE(ex,ey,ez,mu=1):
    '''
        Calculating Electric Field Energy
    '''
    E_E = 0.5*np.mean(ex**2 + ey**2 + ez**2)/mu
    
    return E_E

def Efl(px,py,pz,rm):
    '''
        Calculating Flow Energy from Momentum
    '''
    E_fl = 0.5*np.mean((px**2 + py**2 + pz**2)/rm)
    
    return E_fl

def Eth(Pxx,pyy,pzz,px,py,pz,rm):
    '''
        Calculating Thermal Energy from Pressure Tensor
    '''
    E_th = 0.5*np.mean(Pxx+Pyy+Pzz) - Efl(px,py,pz,rm)
    
    return E_th


    
    
