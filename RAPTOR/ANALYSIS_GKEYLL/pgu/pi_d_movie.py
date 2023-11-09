import numpy as np
import matplotlib.pyplot as plt
import raptor 

N = 250
mr = [25,250,500,1000,1836]
jz = np.zeros((len(mr),N))
c = ['b','g','r','c','m','y','k']
l = [(0, (5, 10)),':','-.','--','-']  #[(0, (1, 10)),(0, (3, 10, 1, 10)),(0, (5, 10)),':','-.','--','-']
jz = np.zeros((5,250,550,550,1))
L = 8*np.pi
ex, ey = L/2, L/2
alpha = np.zeros(5)


t=np.linspace(0,6,250)

for i in range(len(mr)):
    print(mr[i])
    m = mr[i]
    for j in range(0,N):
        print(j)
        jxl, jyl, jzl = raptor.j_xyz(m,j)
        jz[i][j] = jzl
    alpha[i] = np.max(jz[i])

for q in range(len(mr)):
    m = mr[q]
    print(m)
    for k in range(N):
        print(k)
        plt.imshow(jz[q][k]/alpha[q],cmap='seismic',vmin=-1,vmax=1,extent=[-ex,ex,-ey,ey],origin='lower',aspect=1)
        plt.title('$J_z$'+' for mr '+ str(m) +' at timeslice_{}'.format(k))
        plt.colorbar()
#       plt.show()

        plt.savefig('jz_mr'+ str(m) +'_{0:03d}.png'.format(k))

        plt.clf()









