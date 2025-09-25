'''
  _____     .       .    _____  _______  _____   _____
 |//// \   /'(  _  )'\  |//// \|//////_|///// \ |//// \
 |/|__) | / . \/v\/ . \ |/|__) |  |/|  /////   \|/|__) |
 |//_  / /  _)_`-'_(_  \|//___/   |/| |///(_)   |//_  /
 |/| \ \/.-~   ).(   ~-.\/|       |/|  \//     /|/| \ \
 |_|  \/'     /\_/\     `\|       |_|   \_____/ |_|  \_\
              "' '"                                    

Rapid Analyisis for Plasma Turbulence Or Raptor     


PYTHON CODE DEFINING FUNCTIONS TO BE USED IN ANALYSIS OF TURBULENT PLASMA SIMULATIONS RAN IN GKEYLL

'''

import numpy as np
import sys
import os
import scipy.ndimage as ndimage
import eyrie
from eyrie import pdx
from eyrie import pdy
from eyrie import pdz
from eyrie import germano
from eyrie import favre

#Phi_alpha^uT = (germano(P) dot nabla ) dot favre(u)
def Phi_uT(ux,uy,uz,Pxx,Pxy,Pxz,Pyx,Pyy,Pyz,Pzx,Pzy,Pzz,rm,l):
    #calculating u_i u_j
    uxx, uxy, uxz = ux*ux, ux*uy, ux*uz
    uyx, uyy, uyz = uxy  , uy*uy, uy*uz
    uzx, uzy, uzz = uxz  , uyz  , uz*uz
    #calculating thm pressure tensor = P - n*m*u_i*u_j
    pxx, pxy, pxz = Pxx-(rm*uxx), Pxy-(rm*uxy), Pxz-(rm*uxz)
    pyx, pyy, pyz = pxy         , Pyy-(rm*uyy), Pyz-(rm*uyz)
    pzx, pzy, pzz = pxz         , pyz         , Pzz-(rm*uzz)    
    #defining each component of the symetric pressure tensor
    gpxx, gpxy, gpxz = germano(pxx,l), germano(pxy,l), germano(pxz,l) 
    gpyx, gpyy, gpyz = gpxy,           germano(pyy,l), germano(pyz,l)
    gpzx, gpzy, gpzz = gpxz          , gpyz          , germano(pzz,l)  
    #favre filtered u
    fux, fuy, fuz = favre(ux,rm,l), favre(uy,rm,l), favre(uz,rm,l)
    #calculating Phi_alpha^uT
    pgu = -1*(  (gpxx*pdx(fux,dx)) + (gpxy*pdx(fuy,dx)) + (gpxz*pdx(fuz,dx)) \
              + (gpyx*pdy(fux,dy)) + (gpyy*pdy(fuy,dy)) + (gpyz*pdy(fuz,dy)) ) 
    #returning hi_alpha^uT          
    return pgu

#favre filtered tau^u_alpha = favre(u_i u_j) - favre(u_i)favre(u_j) 
def tau_u(ux,uy,uz,rm,l):   	
    #favre filtered u components
    fux, fuy, fuz = favre(ux,rm,l), favre(uy,rm,l), favre(uz,rm,l) 
    #calculating u_i u_j
    uxx, uxy, uxz = ux*ux, ux*uy, ux*uz
    uyx, uyy, uyz = uxy  , uy*uy, uy*uz   
    uzx, uzy, uzz = uxz  , uyz  , uz*uz 
    #filtered u_i u_j: favre(u_i u_j)
    fuxx, fuxy, fuxz = favre(uxx,rm,l), favre(uxy,rm,l), favre(uxz,rm,l)
    fuyx, fuyy, fuyz = fuxy           , favre(uyy,rm,l), favre(uyz,rm,l)
    fuzx, fuzy, fuzz = fuxz           , fuyz           , favre(uzz,rm,l)
    #filtered u_i filtered u_j: favre(u_i)favre(u_j)
    fxfx, fxfy, fxfz = fux*fux, fux*fuy, fux*fuz
    fyfx, fyfy, fyfz = fxfy   , fuy*fuy, fuy*fuz
    fzfx, fzfy, fzfz = fxfz   , fyfz   , fuz*fuz
    #calculating tau_u: favre(u_i u_j) - favre(u_i)favre(u_j)
    txx, txy, txz = fuxx - fxfx, fuxy - fxfy, fuxz - fxfz
    tyx, tyy, tyz = txy        , fuyy - fyfy, fuyz - fyfz
    tzx, tzy, tzz = txz        , tyz        , fuzz - fzfz
    #return tau_u components
    return txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz
    
#calculating favre filtered tau^e = favre(E) - germano(E)
def tau_e(ex,ey,ez,rm,l):
    #calculating favre and germano filtered E
    fex, fey, fez = favre(ex,rm,l), favre(ey,rm,l), favre(ez,rm,l)
    gex, gey, gez = germano(ex,l), germano(ey,l), germano(ez,l)
    #calculating favre filtered tau^e = favre(E) - germano(E)
    tex, tey, tez = fex - gex, fey - gey, fez - gez
    #returning tau^e components
    return tex, tey, tez
    
#calculating favre filtered tau^b = favre(u X B) - favre(u) X favre(B)
def tau_b(bx,by,bz,ux,uy,uz,rm,l):
    #filtering u and B
    fbx, fby, fbz = favre(bx,rm,l), favre(by,rm,l), favre(bz,rm,l)
    fux, fuy, fuz = favre(ux,rm,l), favre(uy,rm,l), favre(uz,rm,l)
    #calculating u X B
    uXBx, uXBy, uXBz = (uy*bz) - (uz*by), -1*((ux*bz) - (uz*bx)), (ux*by) - (uy*bx)
    #calculating favre filtered u X B: favre(u X B)
    fCx, fCy, fCz = favre(uXBx,rm,l), favre(uXBy,rm,l), favre(uXBz,rm,l)
    #calculating favre(u) X favre(B)
    Cfx, Cfy, Cfz = (fuy*fbz) - (fuz*fby), -1*((fux*fbz) - (fuz*fbx)), (fux*fby) - (fuy*fbx)
    #calculating favre filtered tau^b = favre(u X B) - favre(u) X favre(B)
    tbx, tby, tbz = fCx - Cfx, fCy - Cfy, fCz - Cfz
    #return tau^b components
    return tbx, tby, tbz
    
#calculating Pi_alpha^bb = -q germano(n) tau_e dot favre(u)
def Pi_bb(ex,ey,ez,ux,uy,uz,rm,n,q,l):
    #calculating germano(n) and favre(u)
    gn = germano(n,l)
    fux, fuy, fuz = favre(ux,rm,l), favre(uy,rm,l), favre(uz,rm,l)
    #caluculating -q germano(n) tau_e
    tex, tey, tez = tau_e(ex,ey,ez,rm,l)
    qntx, qnty, qntz = -1*q*gn*tex, -1*q*gn*tey, -1*q*gn*tez
    #calculating Pi_alpha^bb = -q germano(n) tau_e dot favre(u)
    Pib = qntx*fux + qnty*fuy + qntz*fuz
    #return Pi_alpha^bb
    return Pib
   
#Pi_alpha^uu = -(german(rho) tau_b dot nabla) dot favre(u) - q/c germano(n) tau_b dot favre(u)
def Pi_uu(bx,by,bz,ux,uy,uz,rm,n,dx,dy,dz,q,l,c=1):
    #defining germano(rho), tau_u, favre(u), germano(n), tau_b
    gr, gn = germano(rm,l), germano(n,l)
    fux, fuy, fuz = favre(ux,rm,l), favre(uy,rm,l), favre(uz,rm,l)
    tuxx, tuxy, tuxz, tuyx, tuyy, tuyz, tuzx, tuzy, tuzz = tau_u(ux,uy,uz,rm,l)
    tbx, tby, tbz = tau_b(bx,by,bz,ux,uy,uz,rm,l)
    #calculating -(german(rho) tau_b dot nabla) dot favre(u)
    tgu = -1*gr*(  (tuxx*pdx(fux,dx)) + (tuxy*pdx(fuy,dx)) + (tuxz*pdx(fuz,dx)) \
                 + (tuyx*pdy(fux,dy)) + (tuyy*pdy(fuy,dy)) + (tuyz*pdy(fuz,dy)) )
    #calculating q/c germano(n) tau_b dot favre(u)
    tdu = (q/(c*gn))*(tbx*fux + tby*fuy + tbz*fuz)
    #calculating Pi_alpha^uu 
    Pi_uu = tgu - tdu
    #returning Pi_alpha_uu
    return Pi_uu

#Favre filtered fluid flow energy: favre(E)^f = 1/2*germano(rho)*favre(u)^2 
def filtered_Efl(ux,uy,uz,rm,l):
    #calculating germano(rho) and favre(u)
    gr = germano(rm,l)
    fux, fuy, fuz = favre(ux,rm,l), favre(uy,rm,l), favre(uz,rm,l)
    Efl = 0.5*(gr*(fux**2 + fuy**2 + fuz**2))
    return Efl

#Germano filtered electromagnetic energy: germano(E)^m = (germano(B)^2+germano(E)^2)/(8*pi)
def filtered_EM(bx,by,bz,ex,ey,ez,l):
    #calculating germano(B) and germano(E)
    gbx, gby, gbz = germano(bx,l), germano(by,l), germano(bz,l)
    gex, gey, gez = germano(ex,l), germano(ey,l), germano(ez,l)
    EM = (((gbx**2 + gby**2 + gbz**2)) + ((gex**2 + gey**2 + gez**2)))/(8*np.pi)
    #return germano(E^m)
    return EM

#divergence of heatflux dmqijm; \partial_{m} Q_{ijm} = v_{th} abs(k_0) (P_{ij} - p\delta_{ij})
def Qij(ux,uy,uz,Pxx,Pxy,Pxz,Pyx,Pyy,Pyz,Pzx,Pzy,Pzz,rm,n,m,q,k0,kB=1):
    #calculating u_i u_j
    uxx, uxy, uxz = ux*ux, ux*uy, ux*uz
    uyx, uyy, uyz = uxy  , uy*uy, uy*uz
    uzx, uzy, uzz = uxz  , uyz  , uz*uz
    #calculating thm pressure tensor
    pxx, pxy, pxz = Pxx-(rm*uxx), Pxy-(rm*uxy), Pxz-(rm*uxz)
    pyx, pyy, pyz = pxy         , Pyy-(rm*uyy), Pyz-(rm*uyz)
    pzx, pzy, pzz = pxz         , pyz         , Pzz-(rm*uzz)
    #calculating scalar pressure
    spe = (1/3)*(pexx+peyy+pezz)
    #calculating temperatrue and thermal velocity
    T = spe/(n*kB)
    vth = np.sqrt(T/m)
    #calculating the deviatoric pressure tensor
    Dxx, Dxy, Dxz = pexx - spe, pexy, pexz
    Dyx, Dyy, Dyz = peyx, peyy - spe, peyz
    Dzx, Dzy, Dzz = pezx, pezy, pezz - spe
    #calculating divergence of heatflux; partial Qijm / partial m
    Qexx, Qexy, Qexz = vth*Dxx*abs(k0), vth*Dxy*abs(k0), vth*Dxz*abs(k0)
    Qeyx, Qeyy, Qeyz = vth*Dyx*abs(k0), vth*Dyy*abs(k0), vth*Dyz*abs(k0)
    Qezx, Qezy, Qezz = vth*Dzx*abs(k0), vth*Dzy*abs(k0), vth*Dzz*abs(k0)
    #returning divergence of heatflux
    return Qexx, Qexy, Qexz, Qeyx, Qeyy, Qeyz, Qezx, Qezy, Qezz, Qixx, Qixy, Qixz, Qiyx, Qiyy, Qiyz, Qizx, Qizy, Qizz

