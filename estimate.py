import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# d=1
y = np.array([0.15015015015015015, 0.41608040201005025, 0.9867886178861789, 2.224469160768453, 4.599394550958627, 9.13502538071066, 17.421159715157682, 29.672782874617738, 45.30541368743616, 64.5316973415133, 85.76428571428572])
x = np.linspace(0.1, 0.3, num=len(y))

for i in range(len(x)):
    y[i] = y[i] / (2500.0 - 2500.0 * x[i]) * 1000.0

plt.plot(x, y, 'b-', label='data')

def func(x, a, b, c):
    return a * np.exp(b * x) + c

popt, pcov = curve_fit(func, x, y, maxfev=5000)
plt.plot(x, func(x, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.legend()
plt.show()

print(popt)
