import numpy as np
from matplotlib import pyplot as plt

data=np.asarray(np.loadtxt("echo2D.out"))
data1=np.asarray(np.loadtxt("echo2D_withJ.out"))
data=np.transpose(data)
data1=np.transpose(data1)

plt.plot(data[0], data[1])
plt.plot(data[0], data1[1])

plt.legend()
plt.show()
