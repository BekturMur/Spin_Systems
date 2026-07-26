import numpy as np
import math as m
import scipy, scipy.integrate, scipy.special, scipy.interpolate
from pylab import *
import sys
import zipfile
import glob
import os
import cmath as CM
import time
import string
import bisect
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.ticker import FormatStrFormatter

import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('echo2Dbezsqrt.out','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(row[0])
        y.append(row[1])

plt.plot(x,y, label='Loaded from file!')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interesting Graph\nCheck it out')
plt.legend()
plt.show()
