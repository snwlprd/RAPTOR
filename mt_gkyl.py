#!/usr/bin/env python
###########################################################################
##
##
###########################################################################
'''
                     _          
                    //\
                   ////\
                  //////\                   /\
                 //\/\/\/\      /\         ///\
     /\         /'(  _  )'\    ///\       /////\
    ///\   /\  / . \/v\/ . \  /\/\/\     /\/\/\/\
   /////\ ///\/  _)_`-'_(_  \/      \   /        \
  /\/\/\/\   /.-~   ).(   ~-.\       \ /          \
 /        \ /'     /\_/\     `\       /__  _  _   _\_     
            _______"' '" ______\     |//_||/| \/_//|/|                 
          ////////////|//////////   /|/|_ |/// |/| |/|      
         /////\///\///|   |//|   \ / |_//||/_\ |/| |/|_      
        /////  \/  |//|   |//|    \_   |_||_|_\|_| | __|         
       /////       |//|   |//|   |//\                   \
      / __/        |__|   |__|   | __\                   \

'''

import numpy as np
import postgkyl as pg
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from os.path import basename, realpath, exists
for dd in ['WorkSpace', 'WorkSpace/python', 'AJGAR/TurbAn']:
    if os.path.exists(os.environ['HOME']+'/'+dd):
        sys.path.insert(0,os.environ['HOME']+'/'+dd+'/')
for dd in ['sims/MR_test_closure','sims/MR_test_workload', 'sims/MR_test','sims/MR_sims', 'analysis/MR_analysis', 'analysis/init']:
    if os.path.exists('/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/'):
        sys.path.insert(0,'/nfs/scratch/edyveaja/GYKL_RUNS'+'/'+dd+'/')
import scipy.ndimage as ndimage
import lpg


def load_data(d, j, claws):
    '''
        Loading in the data from the simulations
    '''
    if isinstance(d, int):
        data = lpg.lpg(str(d))
        data.vars2load(claws)
        data.loadslice(j)
    if isinstance(d, dict):
        data = {}
        for i in range(len(dd)):
            dl = lpg.lpg(str(dd[i]))
            dl.vars2load(claws)
            dl.loadslice(j)
            data[i] = {str(dd[i]):dl}
    if isinstance(d, str):
        data = lpg.lpg(d)
        data.vars2load(claws)
        data.loadslice(j)
    else:
        raise ValueError('Bad value %s (must be int, str or dict)' % mr)
    return data

def B(d,j):
    '''
        Loading magnetic field components bx,by,bz
    '''
    data = load_data(d,j,('bx','by','bz'))
    B = np.array([data.bx,data.by,data.bz])
    return B

def E(d,j):
    '''
        Loading electric field components ex,ey,ez
    '''
    data = load_data(d,j,('ex','ey','ez'))
    E = np.array([data.ex,data.ey,data.ez])
    return E

def p_e(d,j):
    '''
        Loading electron momentum components px,py,pz; p = mnu
    '''
    data = load_data(d,j,('jex','jey','jez'))
    p = np.array([data.jex,data.jey,data.jez])
    return p

def p_i(d,j):
    '''
        Loading ion momentum components px,py,pz; p = mnu
    '''
    data = load_data(d,j,('jix','jiy','jiz'))
    p = np.array([data.jix,data.jiy,data.jiz])
    return p

def p(d,j):
    '''
        Loading electron and ion momentum components px,py,pz; p = mnu
    '''
    data = load_data(d,j,('jex','jey','jez','jix','jiy','jiz'))
    p_e, p_i = np.array([data.jex,data.jey,data.jez]), np.array([data.jix,data.jiy,data.jiz])
    return p_e, p_i

def P_e(d,j):
    '''
        Loading electron (symetric) pressure tensor components Pxx,Pxy,Pxz,Pyy,Pyz,Pzz
    '''
    data = load_data(d,j,('pexx','pexy','pexz','peyy','peyz','pezz'))
    P_ij = np.array([data.pexx,data.pexy,data.pexz,data.peyy,data.peyz,data.pezz])
    return P_ij

def P_i(d,j):
    '''
        Loading ion (symetric) pressure tensor components Pxx,Pxy,Pxz,Pyy,Pyz,Pzz
    '''
    data = load_data(d,j,('pixx','pixy','pixz','piyy','piyz','pizz'))
    P_ij = np.array([data.pixx,data.pixy,data.pixz,data.piyy,data.piyz,data.pizz])
    return P_ij

def P(d,j):
    '''
        Loading electron and ion (symetric) pressure tensor components 
    '''
    data = load_data(d,j,('pexx','pexy','pexz','peyy','peyz','pezz','pixx','pixy','pixz','piyy','piyz','pizz'))
    Pe, Pi = np.array([data.pexx,data.pexy,data.pexz,data.peyy,data.peyz,data.pezz]), np.array([data.pixx,data.pixy,data.pixz,data.piyy,data.piyz,data.pizz])
    return Pe, Pi

def Diag_P_e(d,j):
    '''
        Loading diagonal elctron pressure tensor components Pxx, Pyy, Pzz
    '''
    data = load_data(d,j,('pexx','peyy','pezz'))
    P_ii = np.array([data.pexx,data.peyy,data.pezz])
    return P_ii


def Diag_P_i(d,j):
    '''
        Loading diagonal ion pressure tensor components Pxx, Pyy, Pzz
    '''
    data = load_data(d,j,('pixx','piyy','pizz'))
    P_ii = np.array([data.pixx, data.piyy, data.pizz])
    return P_ii

def Diag_P(d,j):
    '''
        Loading electron and ion pressure tensor components Pxx, Pyy, Pzz
    '''
    data = load_data(d,j,('pexx','peyy','pezz','pixx','piyy','pizz'))
    Pe,Pi = np.array([data.pexx,data.peyy,data.pezz]), np.array([data.pixx,data.piyy,data.pizz])
    return Pe, Pi





      
