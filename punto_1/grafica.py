import numpy as np
import matplotlib.pyplot as plt
import corner

inp=np.genfromtxt('M.txt')
M=int(inp[0])
scale=inp[1]
g=np.genfromtxt('datos.dat', delimiter = '\n')

xs=g[0:M]
ys=g[M:2*M]

params=np.zeros((M,2))
params[:,0]=xs*scale
params[:,1]=ys*scale

fig=plt.figure()
ax=plt.axes()
corner.corner(params, title_kwargs={"fontsize": 12}, show_titles=True, labels=[r'$x$', r'$y$'])
plt.savefig('params.pdf', format='pdf')
plt.close()
