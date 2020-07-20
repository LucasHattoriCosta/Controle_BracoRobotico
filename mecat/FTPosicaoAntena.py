import numpy as np
import scipy as sc
import control as co_general
import matplotlib.pyplot as plt
import control.matlab as co
# Fecha todas as janelas
plt.close('all')
# Parametros da funcao de transferencia
K1 = 100;
Km = 2.083;
Kg = 0.1;
a  = 100;
am = 1.71;
#
# Funcao de transferencia da posicao angular do sistema
s = co.tf('s');
Gtheta = (K1*Km*Kg)/(s*(s+a)*(s+am))
#
# Calculo dos polos da funcao de transferencia Gm
print('-------------')
print('POLOS DA FUNCAO DE TRANSFERENCIA')
print(co.pole(Gtheta))
print('-------------')
print('FT DA PLANTA Gtheta(s) = ')
print(Gtheta)

#
# Expansao em fracoes parciais
# caso nao haja multiplos polos:
#       B(s)       R(1)       R(2)             R(n)
#       ----  =  -------- + -------- + ... + -------- + K(s)
#       A(s)     s - P(1)   s - P(2)         s - P(n)
#
print('-------------')
print('EXPANSAO EM FRACOES PARCIAIS')
B = K1*Km*Kg;
A = [1, a+am, a*am, 0];
[R,P,K]=sc.signal.residue(B,A)
print('R =',R)
print('P =',P)
print('K =',K)
# Plot dos polos e zeros PZMAP
co.pzmap(Gtheta,grid=True)
# Resposta Transitoria
S=co.stepinfo(Gtheta)
print('-------------')
print('CARACTERISTICAS DA RESPOSTA TRANSITORIA DO SISTEMA')
print('tempo de subida tr = ','%.2f' % S['RiseTime'],'seg')
print('tempo de acomodacao ts = ','%.2f' % S['SettlingTime'],'seg')
print('maximo sobresinal Mp = ',S['Overshoot'])
print('valor de pico thetaomax = ','%.2f' % S['Peak'])
print('instante de pico tp = ','%.2f' % S['PeakTime'],'seg')
print('valor de regime estacionario thetaoss = ','%.2f' % S['SteadyStateValue'])
# Resposta a degrau da planta
plt.figure(2)
thetao, t = co.step(Gtheta)
plt.plot(t,thetao)
plt.title('Resposta a degrau da planta')
plt.xlabel('tempo (s)')
plt.ylabel('posicao angular')
plt.grid()

