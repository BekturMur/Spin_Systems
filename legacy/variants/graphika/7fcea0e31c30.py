import numpy as np
from matplotlib import pyplot as plt

data=np.asarray(np.loadtxt("echo2D.out"))
data=np.transpose(data)
x=np.linspace(0,1,1000)
y=10*np.exp(-x**(2/3))

plt.plot(data[0], data[1],'.')
plt.plot(x,y)


plt.show()
