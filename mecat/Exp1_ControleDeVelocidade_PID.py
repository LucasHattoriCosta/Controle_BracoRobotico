import numpy as np
import control as co_general
import matplotlib.pyplot as plt
import control.matlab as co
# Fecha todas as janelas
plt.close('all')
# Controle de velocidade angular do motor eletrico CC
# Testes com controle PID
#
# Definicao dos valores dos parametros do sistema
#
# Constante do tacometro
Ktac = 0.48;
# Parametros da funcao de transferencia Gw(s)
K1 = 100;
Km = 2.083;
Kg = 0.1;
a  = 100;
am = 1.71;
#
# Funcao de transferencia da velocidade angular do sistema
s = co.tf('s');
Gw = (K1*Km*Kg)/((s+a)*(s+am))
print('\n\n-------------------------')
print('RELATORIO')
print('Funcao de transferencia do sistema')
print('Gw = ',Gw)
#
# Controlador PID H1(s)
#                      1        Td*s
#   H1(s) = Kp*( 1 + ----- + ------------ )    
#                     Ti*s     (Td/N)*s+1
# 
#   A parte derivativa possui um Filtro de 1a. ordem
#   com polo em s = (-N/Td)
#
# Qto Maior Ti menor o seu efeito 
# Qto maior N menor o efeito do Filtro de 1a. Ordem
#

# Controlador Proporcional
#H1 = Kp
# Controlador PI
#H1 = Kp*(1+ 1/(Ti*s))
# Controlador PD sem filtro
#H1 = Kp*(1+Td*s)
# Controlador PD com filtro
#H1 = Kp*(1+Td*s/(Td*s/N+1))

# Definicao dos controladores
# Escolha dos parametros
#
# Controlador 1
Kp1 = 1;
Ti1 = 1;
Td1 = 1;
N1  = 1;
# Escolher um dos controladores
# Controlador Proporcional
H1 = Kp1
# Controlador PI
#H1 = Kp1*(1+ 1/(Ti1*s))
# Controlador PD sem filtro
#H1 = Kp1*(1+Td1*s)
# Controlador PD com filtro
#H1 = Kp1*(1+Td1*s/(Td1*s/N1+1))
Hw1 = Ktac*H1
# Controlador 2
Kp2 = 1;
Ti2 = 1;
Td2 = 1;
N2  = 1;
# Escolher um dos controladores
# Controlador Proporcional
H2 = Kp2
# Controlador PI
#H2 = Kp2*(1+ 1/(Ti2*s))
# Controlador PD sem filtro
#H2 = Kp2*(1+Td2*s)
# Controlador PD com filtro
#H2 = Kp2*(1+Td2*s/(Td2*s/N2+1))
Hw2 = Ktac*H2
# Controlador 3
Kp3 = 1;
Ti3 = 1;
Td3 = 1;
N3  = 1;
# Escolher um dos controladores
# Controlador Proporcional
H3 = Kp3
# Controlador PI
#H3 = Kp3*(1+ 1/(Ti3*s))
# Controlador PD sem filtro
#H3 = Kp3*(1+Td3*s)
# Controlador PD com filtro
#H3 = Kp3*(1+Td3*s/(Td3*s/N3+1))
Hw3 = Ktac*H3
#
# Definicao da malha aberta
# Sao definidos 3 sistemas distintos
#
GHw1 = Hw1*Gw
GHw2 = Hw2*Gw
GHw3 = Hw3*Gw
print('\nFUNCAO DE TRANSFERENCIA DE MALHA ABERTA')
print('GHw1 = ',GHw1)
print('GHw2 = ',GHw2)
print('GHw3 = ',GHw3)
#
# Polos e zeros de malha aberta
#
print('\nPOLOS E ZEROS DE MALHA ABERTA')
print('-------------')
print('POLOS E ZEROS - GHw1')
print('Polos = ',co.pole(GHw1))
print('Zeros = ',co.zero(GHw1))
print('-------------')
print('POLOS E ZEROS - GHw2')
print('Polos = ',co.pole(GHw2))
print('Zeros = ',co.zero(GHw2))
print('-------------')
print('POLOS E ZEROS - GHw3')
print('Polos = ',co.pole(GHw3))
print('Zeros = ',co.zero(GHw3))
#
# definicao da malha fechada
# realimentacao unitaria
#
# Funcao de transferencia em malha fechada
# pode ser definida em Matlab utilizando-se
# O comando feedback(S1,S2) 
# Onde o sistema S1 se encontra na malha direta
# e S2 se encontra na malha de realimentacao
# No caso desse sistema a malha direta
# e' G(s)H(s) e malha de realimentacao e' 
# unitaria
#
# R(s)  E(s)|------|  |------|  U(s)
#---->(+)---| H(s) |--| G(s) |------->
#    _ ^    |------|  |------|    |
#      |---------------------------                       
#      
cloop1 = co.feedback(GHw1,1)
cloop2 = co.feedback(GHw2,1)
cloop3 = co.feedback(GHw3,1)
print('\nFUNCAO DE TRANSFERENCIA DE MALHA FECHADA')
print('GHw1/(1+GHw1) = ',cloop1)
print('GHw2/(1+GHw2) = ',cloop2)
print('GHw3/(1+GHw3) = ',cloop3)
#
# Polos e zeros de malha fechada
# comando pole() obtem os polos do sistema
# e zero() obtem os zeros do sistema
#
print('\nPOLOS E ZEROS DE MALHA FECHADA')
print('-------------')
print('Polos e zeros cloop1')
print('Polos = ',co.pole(cloop1))
print('Zeros = ',co.zero(cloop1))
print('COEF. DE AMORTECIMENTO E FREQ. NATURAL')
print('_____Polos____________zeta_______omegan')
co.damp(cloop1)
print('-------------')
print('Polos e zeros cloop2')
print('Polos = ',co.pole(cloop2))
print('Zeros = ',co.zero(cloop2))
print('COEF. DE AMORTECIMENTO E FREQ. NATURAL')
print('_____Polos____________zeta_______omegan')
co.damp(cloop2)
print('-------------')
print('Polos e zeros cloop3')
print('Polos = ',co.pole(cloop3))
print('Zeros = ',co.zero(cloop3))
print('COEF. DE AMORTECIMENTO E FREQ. NATURAL')
print('_____Polos____________zeta_______omegan')
co.damp(cloop3)
#
# Grafico da resposta a degrau unitario
#
plt.figure(1)
# Definicao do vetor de tempo t
t=[x*0.2 for x in range(0,50)]
omegao1, t = co.step(cloop1,t)
omegao2, t = co.step(cloop2,t)
omegao3, t = co.step(cloop3,t)
plt.plot(t,omegao1,t,omegao2,t,omegao3)
plt.title('Resposta a degrau do sistema em malha fechada')
plt.xlabel('tempo (s)')
plt.ylabel('velocidade Angular')
plt.grid()
# Calcula as caracteristicas da resposta transitoria
#  stepinfo(sys, T=None, SettlingTimeThreshold=0.02, RiseTimeLimits=(0.1,0.9))
#  S: a dictionary containing:
#        RiseTime: Time from 10% to 90% of the steady-state value.
#        SettlingTime: Time to enter inside a default error of 2%
#        SettlingMin: Minimum value after RiseTime
#        SettlingMax: Maximum value after RiseTime
#        Overshoot: Percentage of the Peak relative to steady value
#        Undershoot: Percentage of undershoot
#        Peak: Absolute peak value
#        PeakTime: time of the Peak
#        SteadyStateValue: Steady-state value
S1 = co.stepinfo(cloop1)
print('-------------')
print('CARACTERISTICAS DA RESPOSTA TRANSITORIA DO SISTEMA cloop1')
print('tempo de subida tr = ','%.2f' % S1['RiseTime'],'seg')
print('tempo de acomodacao ts = ','%.2f' % S1['SettlingTime'],'seg')
print('maximo sobresinal Mp = ',S1['Overshoot'])
print('valor de pico thetaomax = ','%.2f' % S1['Peak'])
print('instante de pico tp = ','%.2f' % S1['PeakTime'],'seg')
print('valor de regime estacionario omegaoss = ','%.2f' % S1['SteadyStateValue'])
S2 = co.stepinfo(cloop2)
print('-------------')
print('CARACTERISTICAS DA RESPOSTA TRANSITORIA DO SISTEMA cloop2')
print('tempo de subida tr = ','%.2f' % S2['RiseTime'],'seg')
print('tempo de acomodacao ts = ','%.2f' % S2['SettlingTime'],'seg')
print('maximo sobresinal Mp = ',S2['Overshoot'])
print('valor de pico thetaomax = ','%.2f' % S2['Peak'])
print('instante de pico tp = ','%.2f' % S2['PeakTime'],'seg')
print('valor de regime estacionario omegaoss = ','%.2f' % S2['SteadyStateValue'])
S3 = co.stepinfo(cloop3)
print('-------------')
print('CARACTERISTICAS DA RESPOSTA TRANSITORIA DO SISTEMA cloop3')
print('tempo de subida tr = ','%.2f' % S3['RiseTime'],'seg')
print('tempo de acomodacao ts = ','%.2f' % S3['SettlingTime'],'seg')
print('maximo sobresinal Mp = ',S3['Overshoot'])
print('valor de pico thetaomax = ','%.2f' % S3['Peak'])
print('instante de pico tp = ','%.2f' % S3['PeakTime'],'seg')
print('valor de regime estacionario omegaoss = ','%.2f' % S3['SteadyStateValue'])
#
# Funcao de transferencia para calculo do esforco de controle u(t)
# o sinal u(t) pode ser calculado definindo-se
# um sistema de controle em malha fechada onde H(s)
# esta na malha direta e G(s) na malha de realimentacao
#
# R(s)  E(s)|------|        U(s)
#---->(+)---| H(s) |------------>
#    _ ^    |------|    |
#      |                |
#      |    |------|    |
#      |----| G(s) |<----
#           |------|
#      
esforco1 = co.feedback(Hw1,Gw)
esforco2 = co.feedback(Hw2,Gw)
esforco3 = co.feedback(Hw3,Gw)
vp1, t = co.step(esforco1,t)
vp2, t = co.step(esforco2,t)
vp3, t = co.step(esforco3,t)
plt.figure(2)
plt.plot(t,vp1,t,vp2,t,vp3)
plt.title('Esforco de controle')
plt.xlabel('tempo (s)')
plt.ylabel('vp(t)')
plt.grid()
#
# Caracteristicas do esforco de controle
#
print('\nESFORCO DE CONTROLE')
SE1 = co.stepinfo(esforco1)
print('-------------')
print('CARACTERISTICAS DO ESFORCO DE CONTROLE cloop1')
print('valor de pico vpmax = ','%.2f' % SE1['Peak'])
print('instante de pico tp = ','%.2f' % SE1['PeakTime'],'seg')
print('valor de regime estacionario vpss = ','%.2f' % SE1['SteadyStateValue'])
SE2 = co.stepinfo(esforco2)
print('-------------')
print('CARACTERISTICAS DO ESFORCO DE CONTROLE cloop2')
print('valor de pico vpmax = ','%.2f' % SE2['Peak'])
print('instante de pico tp = ','%.2f' % SE2['PeakTime'],'seg')
print('valor de regime estacionario vpss = ','%.2f' % SE2['SteadyStateValue'])
SE3 = co.stepinfo(esforco3)
print('-------------')
print('CARACTERISTICAS DO ESFORCO DE CONTROLE cloop3')
print('valor de pico vpmax = ','%.2f' % SE3['Peak'])
print('instante de pico tp = ','%.2f' % SE3['PeakTime'],'seg')
print('valor de regime estacionario vpss = ','%.2f' % SE3['SteadyStateValue'])
