import pickle
import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc_int
import subprocess as sub

# # set up latex
# plt.rc('text', usetex=True)

def DoubleIntegrate(z):
	total_sum = np.sum(z)
	return total_sum

vardict={'rme':[0,'_electron_'],'jex':[1,'_electron_'],'jey':[2,'_electron_'],\
       'jez':[3,'_electron_'],'pexx':[4,'_electron_'],'pexy':[5,'_electron_'],\
       'pexz':[6,'_electron_'],'peyy':[7,'_electron_'],'peyz':[8,'_electron_'],\
       'pezz':[9,'_electron_'],'rmi':[0,'_hydrogen_'],'jix':[1,'_hydrogen_'],\
       'jiy':[2,'_hydrogen_'],'jiz':[3,'_hydrogen_'],'pixx':[4,'_hydrogen_'],\
       'pixy':[5,'_hydrogen_'],'pixz':[6,'_hydrogen_'],'piyy':[7,'_hydrogen_'],\
       'piyz':[8,'_hydrogen_'],'pizz':[9,'_hydrogen_'],'ex':[0,'_field_'],\
       'ey':[1,'_field_'],'ez':[2,'_field_'],'bx':[3,'_field_'],'by':[4,'_field_'],\
       'bz':[5,'_field_']}

def load_data(vardict,basename,N,vrlist):
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

#def's for curl
def pderiv(ar,dx=1,axis=0):
        return (np.roll(ar,1,axis=axis)-\
                np.roll(ar,-1,axis=axis))/(2*dx)
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
