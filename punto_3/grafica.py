gimport numpy as np
import matplotlib.pyplot as plt
import corner

g=np.genfromtxt('datos.dat', delimiter=' ')

alfa=g[:,0]
beta=g[:,1]
gamma=g[:,2]
delta=g[:,3]

params=np.zeros((len(alfa),4))
params[:,0]=alfa
params[:,1]=beta
params[:,2]=gamma
params[:,3]=delta

fig=plt.figure()
ax=plt.axes()
corner.corner(params, title_kwargs={"fontsize": 12}, show_titles=True, labels=[r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\delta$'])
plt.savefig('params.pdf', format='pdf')
plt.close()
