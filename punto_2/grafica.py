import numpy as np
import matplotlib.pyplot as plt
import corner

G=6.67408*10**(-11)

M=np.genfromtxt('M.txt', dtype='int')
g=np.genfromtxt('datos.dat', delimiter='\n')

pends=g[0:M]
cortes=g[M:2*M]
like=g[2*M:3*M]

params=np.zeros((M,2))
params[:,0]=1-pends
params[:,1]=cortes*10-np.log10(G)

fig=plt.figure()
ax=plt.axes()
corner.corner(params, title_kwargs={"fontsize": 12}, show_titles=True, labels=[r'$\alpha$', r'$log_{10}M_{sol}$'])
plt.savefig('params.pdf', format='pdf')
plt.close()
