#Alunos - Reinaldo Coutinho, Marcus Bispo, Felipe Viana e Orlando Nascimento
from numpy import pi, sqrt, abs, linspace
from scipy.integrate import solve_ivp as EDO
from scipy.fft import fft
import matplotlib.pyplot as plt

N = 1000
M = 1
k = 1
k12 = 0.1
x10 = 1
x20 = 0
v10 = v20 = 0
x1 = x10
y0 = (x10, v10, x20, v20)

def mov(t, f, M, k, k12):
	xv = f[1]
	vk = -k*f[0]/M + k12*(f[2]-f[0])/M
	xv2 = f[3]
	vk2 = -k*f[2]/M - k12*(f[2]-f[0])/M
	return (xv, vk, xv2, vk2)

def Sistema (y0, t_span, t_eval, M, k, k12):
	return EDO(mov, y0 = y0, t_span=t_span, t_eval = t_eval, args = (M, k, k12), method='RK45')

#a)
t_eval, Dt = linspace(0, 270, 1000, retstep = True)
t_span = (t_eval[0], t_eval[-1])
solucao = Sistema(y0, t_span, t_eval, M, k, k12)

plt.plot(solucao.t, solucao.y[0], '-r')
plt.show()

plt.plot(solucao.t, solucao.y[2], '-b')
plt.show()
#b)
fourier1 = fft(solucao.y[0])[:N//2]
fourier2 = fft(solucao.y[2])[:N//2]

frequencia = linspace(0, 1/(2*Dt), N//2)

cn1 = 2*abs(fourier1)/N
cn2 = 2*abs(fourier2)/N

v1 = sqrt(k/M)/(2*pi)
v2 = sqrt((2*k12+k)/M)/(2*pi)
print ('analítical v1 = ', str(v1))
print ('analitical v2 = ', str(v2))
plt.plot(frequencia, cn1, '-r')
plt.plot(frequencia, cn2, '-b')
plt.show()