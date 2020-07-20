% Parametros da funcao de transferencia
K1 = 100;
Km = 2.083;
Kg = 0.1;
a  = 100;
am = 1.71;
%
% Funcao de transferencia da posicao angular do sistema
s = tf('s');
Gp = (K1*Km*Kg)/(s*(s+a)*(s+am))
%
% Calculo dos polos da funcao de transferencia Gp
pole(Gp)
%
% Expansao em fracoes parciais
% caso nao haja multiplos polos:
%       B(s)       R(1)       R(2)             R(n)
%       ----  =  -------- + -------- + ... + -------- + K(s)
%       A(s)     s - P(1)   s - P(2)         s - P(n)
%
B=K1*Km*Kg;
A=[1 a+am a*am 0];
[R,P,K]=residue(B,A)
%
% Plot dos polos e zeros de Gm
%
figure(1)
pzmap(Gp)
% Resposta a entrada degrau do sistema
figure(2)
t=0:0.02:10;
step(Gp,t)
grid on
% Caracteristicas da resposta a degrau
stepinfo(Gp)
