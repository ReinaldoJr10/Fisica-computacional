# Grupo : Orlando Nascimento, Marcus Bispo, Reinaldo Coutinho, Felipe Viana
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

q=m=E=B=1

def eletrica_magnetica_edo(t,y):
  dy2_dt2=y[1]
  dvy2_dt2=q*B*y[3]
  dz2_dt2=y[3]
  dvz2_dt2=q*(E-B*(y[1]))
  return [dy2_dt2,dvy2_dt2,dz2_dt2,dvz2_dt2]

inicio = 0
fim    = 50
t_eval = np.linspace (inicio, fim, 1000)
t_span = [min(t_eval ),max(t_eval)] 
# instantes incial e final .

y0 = [0.0,0.0,0.0,0.0]
#Condicoes iniciais: y,vy,z,vz

# integrador de EDO da SciPy .
campo = integrate.solve_ivp (eletrica_magnetica_edo,
t_span = t_span ,
y0 = y0 ,
t_eval = t_eval ,
method='DOP853',
rtol = 1.e-10 ,
atol = 1.e-10)

plt.title ('Gráfico')
plt.plot (campo.y[0], campo.y[2], '-b')
plt.xlabel ('y')
plt.ylabel ('z')
plt.show ()

R = np.sqrt ((campo.y[0] - campo.t)**2 + (campo.y[2] - 1)**2)
plt.plot (campo.t, R, '-b')
plt.axhline (1, ls='--', lw=2)
plt.xlabel ('t(s)')
plt.ylabel ('$R(m)$')
plt.show ()

#a) O Período de repetição é de aproximadamente 6 segundos.

#b) O valor de R encontrado seria como uma espécie de erro relativo nos calculos, o valor varia em valores muitos pequenos ao redor de 1  ...Como y(t)=t-sin(t) e z(t)=1-cos(t) podemos re-escrever as equações como: y(t)-t = -sin(t) e z(t)-1 = -cos(t), substituido esses valores na equação de R = sqrt( (-sin(t))² + (-cos(t))² ) que pode simplificar como R = sqrt( sin²(t) + cos²(t)), usando a identidade trigonometrica sin²(x)+cos²(x)=1, podemos substituir em R = sqrt(1) = 1. Portanto, R(t) = 1.  