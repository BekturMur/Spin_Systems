import numpy as np
from matplotlib import pyplot as plt

data=np.asarray(np.loadtxt("echo2D.out"))
data=np.transpose(data)

plt.plot(data[0], data[1],'.')
plt.xlabel(r'$t$')
plt.ylabel(r'$M_{z}(t)$')


plt.show()
