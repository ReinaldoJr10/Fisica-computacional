#Alunos - Gabriel Ribeiro, Lucas Parente, Pedro Mendes
from tkinter import N
from numpy import pi, sqrt, abs, linspace
from scipy.integrate import solve_ivp as resolveEDO
from scipy.fft import fft
import matplotlib.pyplot as plt

def movimentacao(t, f, M, k, k12):
	xv = f[1]
	vk = -k*f[0]/M + k12*(f[2]-f[0])/M
	xv2 = f[3]
	vk2 = -k*f[2]/M - k12*(f[2]-f[0])/M
	return (xv, vk, xv2, vk2)

def resolveSistema (y0, t_span, t_eval, M, k, k12):
	return resolveEDO(movimentacao, y0 = y0, t_span=t_span, t_eval = t_eval, args = (M, k, k12), method='RK45')

N = 1000

M = 1
k = 1
k12 = 0.1

x10 = 1
x20 = 0
v10 = v20 = 0
x1 = x10
y0 = (x10, v10, x20, v20)

t_eval, Dt = linspace(0, 270, 1000, retstep = True)
t_span = (t_eval[0], t_eval[-1])
solucao = resolveSistema(y0, t_span, t_eval, M, k, k12)

plt.plot(solucao.t, solucao.y[0], '-r')
plt.show()

plt.plot(solucao.t, solucao.y[2], '-b')
plt.show()

fx1 = fft(solucao.y[0])[:N//2]
fx2 = fft(solucao.y[2])[:N//2]

freq = linspace(0, 1/(2*Dt), N//2)

cn1 = 2*abs(fx1)/N
cn2 = 2*abs(fx2)/N

v1 = sqrt(k/M)/(2*pi)
v2 = sqrt((2*k12+k)/M)/(2*pi)
print(v1, v2)

plt.plot(freq, cn1, '-r')
plt.plot(freq, cn2, '-b')
plt.show()