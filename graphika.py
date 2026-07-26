import numpy as np
from matplotlib import pyplot as plt

data=np.asarray(np.loadtxt("CPMGZ_copy18spins64/echo2D.out"))
data=np.transpose(data)

plt.plot(data[0], data[1])

plt.xlabel(r'$t$')
plt.ylabel(r'$\langle Re M_{x} \rangle$')
plt.text(1.5,13,'N=18', fontsize =12)


plt.show()
