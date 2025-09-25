#!/usr/bin/env python
###########################################################################
##
##
###########################################################################

import numpy as np
import postgkyl as pg
from os.path import basename, realpath, exists
   
def imsh(d,varname,**kwargs):
    """
    Expects a dictionary that has 'x','y' variables.
    """
    plt.imshow(d[varname][:,:,0].T,origin='lower',extent=[d['x'].min(),d['x'].max(),d['y'].min(),d['y'].max()],**kwargs)
    plt.colorbar()
def loadComponent(data, comp):
    data_comp = pg.data.select(data, comp=comp)
    return np.asarray(data_comp[1])[:, :, 0]

def getXYCoords(data):
    pos = pg.data.select(data, comp=0)
    return np.asarray(pos[0])

#def for derivative
def pderiv(ar,dx=1,axis=0):
    return (np.roll(ar,1,axis=axis)-\
        np.roll(ar,-1,axis=axis))/(2*dx)

#def's for curl
def pdx(ar,dx):
    return pderiv(ar,dx,0)
def pdy(ar,dy):
    return pderiv(ar,dy,0)
def pdz(ar,dz):
    return pderiv(ar,dz,0)
def curl(ax,ay,az,dx,dy,dz):
    cx = pdy(az,dy)-pdz(ay,dz)
    cy = pdz(ax,dz)-pdx(az,dx)
    cz = pdx(ay,dx)-pdy(ax,dy)
    return cx,cy,cz
def curlx(ay,az,dy,dz):
    return pdy(az,dy)-pdz(ay,dz)
def curly(ax,az,dx,dz):
    return pdz(ax,dz)-pdx(az,dx)
def curlz(ax,ay,dx,dy):
    return pdx(ay,dx)-pdy(ax,dy)

#def's for divergene
def pdx(ar,dx):
    return pderiv(ar,dx,0)
def pdy(ar,dy):
    return pderiv(ar,dy,0)
def pdz(ar,dz):
    return pderiv(ar,dz,0)
def div(ax,ay,az,dx,dy,dz):
    divergence = pdx(ax,dx) + pdy(ay,dy) + pdz(az,dz)
    return divergence

def calc_dep(v2lu,primitives,derived):
   # Find out the dependencies and append them before the variable
   lgcl=[False]*len(v2lu)
   # If there is even a single False, go through the loop
   while not all(lgcl):
   # For all values in the list
      for i in v2lu:
   # Find the index of the item in list
         idx=v2lu.index(i)
   # If it is a primitive item, set the dependency satisfied flag to true
         if i in primitives:
            lgcl[idx]=True
   # If it is a derived quantity, check all its dependencied for existence in the list
         elif i in derived:
   # create empty logical and insertion lists to insert dependencies before item i
            tmplgcl=[]; insrt=[]
            for j in derived[i]:
   # If the dependency of the derived quantity is in list, don't do anything
               if j in v2lu:
                  pass
               else:
   # Else append the insertion list
                  insrt.append(j)
   # If the inserted quantity is a primitive, set the dependency satisfied flag to true, else false
                  if j in primitives:
                     tmplgcl.append(True)
                  else:
                     tmplgcl.append(False)
   # Insert the newly found quantities into the list _before_ i
            v2lu=v2lu[:idx]+insrt+v2lu[idx:]
   # Set the dependency satisfied flag for i to be True
            lgcl[idx]=True
   # Insert the dependency satisfied flags for the newly inserted quantities
            lgcl=lgcl[:idx]+tmplgcl+lgcl[idx:]
   return v2lu
class lpg(object):
   """
     Class to load postgkyl data into an object that works
     with TurbAn codes.
     
     - 2022/02/09: First dirty hack from p3d.py sources by TNP
   """
   def __init__(self,filename):
      import importlib
      init = importlib.import_module(filename)
      for k in init.__dict__.keys():
         self.__dict__[k] = init.__dict__[k]
      
      self.primitives=['bx','by','bz','ex','ey','ez','jix','jiy','jiz','jex','jey'\
      ,'jez','rmi','pixx','piyy','pizz','pixy','pixz','piyz'\
      ,'pexx','pexy','pexz','peyy','peyz','pezz','rme','rho']
      self.derived={\
      'ni':['rmi'], 'ne':['rme'],\
      'jx':['by','bz'],'jy':['bx','bz'],'jz':['bx','by'],\
      'tex':['pexx','ne'],'tey':['peyy','ne'],'tez':['pezz','ne'],\
      'te':['tex','tey','tez'],\
      'tix':['pixx','ni'],'tiy':['piyy','ni'],'tiz':['pizz','ni'],\
      'ti':['tix','tiy','tiz'],\
      'vix':['jix','ni'],'viy':['jiy','ni'],'viz':['jiz','ni'],\
      'vi':['vix','viy','viz'],\
      'vex':['jex','ne'],'vey':['jey','ne'],'vez':['jez','ne'],\
      've':['vex','vey','vez'],\
      'omix':['viy','viz'],'omiy':['viz','vix'],'omiz':['vix','viy'],\
      'omi':['omix','omiy','omiz'],'ensti':['omi'],'pali':['omi'],\
      'omex':['vey','vez'],'omey':['vex','vez'],'omez':['vex','vey'],\
      'ome':['omex','omey','omez'],'enste':['ome'],'pale':['ome'],\
      'omx':['cmy','cmz'],'omy':['cmz','cmx'],'omz':['cmx','cmy'],\
      'om':['omx','omy','omz'], 'enst':['om'], 'pal':['om'],\
      'dui':['vix','viy','viz'],'due':['vex','vey','vez'],\
      'den':['ni','ne'],\
      'cmx':['vix','vex'],'cmy':['viy','vey'],'cmz':['viz','vez'],\
      'zpx':['bx','cmx','den'],'zpy':['by','cmy','den'],'zpz':['bz','cmz','den'],\
      'zmx':['bx','cmx','den'],'zmy':['by','cmy','den'],'zmz':['bz','cmz','den'],\
      'zpzm':['zpx','zpy','zpz','zmx','zmy','zmz']\
      }

      self.allvars=self.primitives+list(self.derived.keys())

####
#### Method to define the variables to load, create variables and
#### open corresponding files.
#### 
   def vars2load(self,v2lu):
      """
         Define the variables to load, define corresponding numpy arrays &
         open the files
      """
      if len(v2lu) == 1:
         if v2lu[0] == 'min':
            self.vars2l=['bx','by','jix','jiy','jz','ni']
         elif v2lu[0] == 'prim':
            self.vars2l=self.primitives
         elif v2lu[0] == 'all':
            self.vars2l=self.allvars
         else:
            # Compute all the dependencies and set variables 2 load 
            self.vars2l=calc_dep(v2lu,self.primitives,self.derived)
      else:
         # Compute all the dependencies and set variables 2 load 
         self.vars2l=calc_dep(v2lu,self.primitives,self.derived)

   def readslice(self,vardict,basename,N,vrlist):
       d = {}
       for v in vrlist:
           # Load in the data
           filename=basename + vardict[v][1] + str(N) + '.bp'
           #print(v,filename)
           data = pg.Data(filename)
           # 
           # pg.data.select outputs a tuple ([x,y],variable)
           d[v] = pg.data.select(data, comp=vardict[v][0])[1]
           
       d['x'], d['y'] = np.asarray(pg.data.select(data,comp=vardict[v][0])[0])
       return d

####
#### Method to load time slices for the loaded variables.
####
   def loadslice(self,it,smth=None):
      """
         Load the variables initialized by self.vars2load()
      """
      for i in self.vars2l:
         if i in self.primitives:
            self.__dict__[i]=self.readslice(self.vardict,self.basename,it,[i])[i]
      for i in self.vars2l:
         if i in self.derived:
           #self.__dict__[i] = self._derivedv(i)
            self._derivedv(i)

      self.mmd={}
      for i in self.vars2l:
         if smth is not None:
            self.__dict__[i]=gf(self.__dict__[i],sigma=smth)
         self.mmd[i]=[self.__dict__[i].min(),self.__dict__[i].max()]
      self.time = it*self.dtmovie

####
#### Method to add attributes to the object
####
   def addattr(self,key,val):
      for i in key:
         print('Adding '+i) 
         self.__dict__[i]=val[key.index(i)]
         if isinstance(val[key.index(i)],np.ndarray):
            self.mmd[i]=[self.__dict__[i].min(),self.__dict__[i].max()]

####
#### Method to compute derived quantities
####
   def _derivedv(self,varname):
      if varname == 'ni'    : self.ni     = self.rmi/self.mi
      if varname == 'ne'    : self.ne     = self.rme/self.me
      if varname == 'jx'    : self.jx     = curlx(self.bz,self.by,self.dy,self.dz)
      if varname == 'jy'    : self.jy     = curly(self.bx,self.bz,self.dz,self.dx)
      if varname == 'jz'    : self.jz     = curlz(self.by,self.bx,self.dx,self.dy)
#  not good - dont uncomment unless solved issue, check mean square codes in /MassRatioSims/...    if varname == 'msj'   : self.msj    = np.mean((self.jx**2+self.jy**2+self.jz**2))    
      if varname == 'tix'   : self.tix    = self.pixx/self.ni
      if varname == 'tiy'   : self.tiy    = self.piyy/self.ni
      if varname == 'tiz'   : self.tiz    = self.pizz/self.ni
      if varname == 'ti'    : self.ti     = (self.tix+self.tiy+self.tiz)/3.
      if varname == 'tex'   : self.tex    = self.pexx/self.ne
      if varname == 'tey'   : self.tey    = self.peyy/self.ne
      if varname == 'tez'   : self.tez    = self.pezz/self.ne
      if varname == 'te'    : self.te     = (self.tex+self.tey+self.tez)/3.
      if varname == 'vix'   : self.vix    = self.jix/self.ni
      if varname == 'viy'   : self.viy    = self.jiy/self.ni
      if varname == 'viz'   : self.viz    = self.jiz/self.ni
      if varname == 'omix'  : self.omix   = curlx(self.viz,self.viy,self.dy,self.dz)
      if varname == 'omiy'  : self.omiy   = curly(self.vix,self.viz,self.dz,self.dx)
      if varname == 'omiz'  : self.omiz   = curlz(self.viy,self.vix,self.dx,self.dy)
      if varname == 'omi'   : pass
      if varname == 'ensti' : self.ensti  = self.omix**2+self.omiy**2+self.omiz**2
      if varname == 'dui'   : self.dui    = div(self.vix,self.viy,self.viz,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'vex'   : self.vex    = -self.jex/self.ne
      if varname == 'vey'   : self.vey    = -self.jey/self.ne
      if varname == 'vez'   : self.vez    = -self.jez/self.ne
      if varname == 'omex'  : self.omex   = curlx(self.vez,self.vey,self.dy,self.dz)
      if varname == 'omey'  : self.omey   = curly(self.vex,self.vez,self.dz,self.dx)
      if varname == 'omez'  : self.omez   = curlz(self.vey,self.vex,self.dx,self.dy)
      if varname == 'ome'   : pass
      if varname == 'enste' : self.enste  = self.omex**2+self.omey**2+self.omez**2
      if varname == 'due'   : self.due    = div(self.vex,self.vey,self.vez,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'cmx'   : self.cmx    = (self.vix+self.me*self.vex)/(1+self.me)
      if varname == 'cmy'   : self.cmy    = (self.viy+self.me*self.vey)/(1+self.me)
      if varname == 'cmz'   : self.cmz    = (self.viz+self.me*self.vez)/(1+self.me)
      if varname == 'omx'   : self.omx    = curlx(self.cmy,self.cmz,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omy'   : self.omy    = curly(self.cmz,self.cmx,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'omz'   : self.omz    = curlz(self.cmx,self.cmy,dx=self.dx,dy=self.dy,dz=self.dz)
      if varname == 'om'    : pass
      if varname == 'enst'  : self.enst   = self.omx**2+self.omy**2+self.omz**2
      if varname == 'den'   : self.den    = self.rmi+self.rme
      if varname == 'zpx'   : self.zpx    = self.bx/np.sqrt(self.den) + self.cmx
      if varname == 'zpy'   : self.zpy    = self.by/np.sqrt(self.den) + self.cmy
      if varname == 'zpz'   : self.zpz    = self.bz/np.sqrt(self.den) + self.cmz
      if varname == 'zmx'   : self.zmx    = self.bx/np.sqrt(self.den) - self.cmx
      if varname == 'zmy'   : self.zmy    = self.by/np.sqrt(self.den) - self.cmy
      if varname == 'zmz'   : self.zmz    = self.bz/np.sqrt(self.den) - self.cmz
      if varname == 'zpzm'  : pass

      if varname == 'pali'  : 
         tmp = af.pcurl(self.omix,self.omiy,self.omiz,dx=self.dx,dy=self.dy,dz=self.dz)
         self.pali   = 0.5*(tmp[0]**2+tmp[1]**2+tmp[2]**2); tmp=None
      if varname == 'pale'  : 
         tmp = af.pcurl(self.omex,self.omey,self.omez,dx=self.dx,dy=self.dy,dz=self.dz)
         self.pale   = 0.5*(tmp[0]**2+tmp[1]**2+tmp[2]**2); tmp=None
      if varname == 'pal'  : 
         tmp = af.pcurl(self.omx,self.omy,self.omz,dx=self.dx,dy=self.dy,dz=self.dz)
         self.pal    = 0.5*(tmp[0]**2+tmp[1]**2+tmp[2]**2); tmp=None
