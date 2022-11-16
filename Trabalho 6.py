import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.fft import fft

N = 1000

def fun (t, y):
    f    = np.empty (y.size)

    f[0] = y[1]
    f[1] = ((-y[5]/y[4])*y[0]) + ((y[5]/y[4])*(y[2]-y[0]))
    f[2] = y[3]
    f[3] = ((-y[5]/y[4])*y[2]) - ((y[5]/y[4])*(y[2]-y[0]))
    return f
M = 1
k = 1
k12 = 0.1
x10 = 1
x20 = 0
v10 = 0
v20 = 0
te, dt  = np.linspace (0., 270., N, retstep=True)
ts  = (te.min(), te.max())
y0  = np.array ([x10, v10, x20, v20, M, k, k12])
res = solve_ivp (fun, t_span=ts, y0=y0, method='RK45', t_eval=te)
plt.plot (res.t, res.y[0], '-r')
plt.plot (res.t, res.y[2], '-b')
plt.show()

ffx = fft (res.y[3])[:N//2]
ffy = fft (res.y[4])[:N//2]
freq = np.linspace(0, 1 /(2 * dt ), N//2)
cn_x = 2 * np.abs (ffx) / N
cn_y = 2 * np.abs (ffy) / N

v1 = (1/(2*np.pi))*(np.sqrt(k/M))
v2 = (1/(2*np.pi))*(np.sqrt(((2*k12)+k)/M))
print ('analítical v1 = ', str(v1))
print ('analitical v2 = ', str(v2))

plt.plot(freq, cn_x, '-r')
plt.plot(freq, cn_y, '-b')
plt.show()