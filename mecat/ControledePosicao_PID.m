%
% Controle de posicao angular do motor eletrico CC
% Testes com controlador PID
%
% Definicao dos valores dos parametros do sistema
%
% Constante do tacometro
Kpot = 0.318;
% Parametros da funcao de transferencia Gw(s)
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
% Controlador PID H(s)
%                     1        Td*s
%   H(s) = Kp*( 1 + ----- + ------------ )    
%                    Ti*s     (Td/N)*s+1
% 
%   A parte derivativa possui um Filtro de 1a. ordem
%   com polo em s = (-N/Td)
%
% Para controlador P:    Ti=Inf, Td=0, N=Inf, Kp > 0
% Para controlador PI:   Td=0, N=Inf, Kp, Ti > 0
% Para controlador PD:   Ti=Inf, Kp, Ti, N > 0
% Para controlador PID:  Kp, Ti, Td, N >0
%
% Qto Maior Ti menor o seu efeito 
% Qto maior N menor o efeito do Filtro de 1a. Ordem
%
% Definicao dos controladores
% Escolha dos parametros
%
% Definicao dos controladores
% Escolha dos parametros
%
% Controlador 1
Kp1 = 10
Ti1 = Inf;
Td1 = 0;
N1  = Inf;
H1  = pidstd(Kp1,Ti1,Td1,N1)
Hp1 = Kpot*H1
% Controlador 2
Kp2 = 20
Ti2 = Inf;
Td2 = 0;
N2  = Inf;
H2  = pidstd(Kp2,Ti2,Td2,N2)
Hp2 = Kpot*H2
% Controlador 3
Kp3 = 50
Ti3 = Inf;
Td3 = 0;
N3  = Inf;
H3  = pidstd(Kp3,Ti3,Td3,N3)
Hp3 = Kpot*H3
%
% Definicao da malha aberta
% Sao definidos 3 sistemas distintos
%
GHp1 = Hp1*Gp
GHp2 = Hp2*Gp
GHp3 = Hp3*Gp
%
% Polos e zeros de maplha aberta
%
display('Polos e zeros - GHp1')
pole(GHp1)
zero(GHp1)
display('Polos e zeros - GHp2')
pole(GHp2)
zero(GHp2)
display('Polos e zeros - GHp3')
pole(GHp3)
zero(GHp3)
%
% Plot com os polos e zeros de malha aberta
%
figure(1);
subplot(1,3,1);
pzplot(GHp1);
title('m. aberta GHp1');
axis equal;
subplot(1,3,2);
pzplot(GHp2);
title('m. aberta GHp2');
axis equal;
subplot(1,3,3);
pzplot(GHp3);
title('m. aberta GHp3');
axis equal;
%
% definicao da malha fechada
% realimentacao unitaria
%
cloop1 = feedback(GHp1,1)
cloop2 = feedback(GHp2,1)
cloop3 = feedback(GHp3,1)
%
% Polos e zeros de malha fechada
%
display('Polos e zeros cloop1')
pole(cloop1)
zero(cloop1)
display('Polos e zeros cloop2')
pole(cloop2)
zero(cloop2)
display('Polos e zeros cloop3')
pole(cloop3)
%
%  Plot com os polos e zeros de malha fechada
%
figure(2);
subplot(1,3,1)
pzplot(cloop1);     
title('m. fechada CL1');
axis equal;
subplot(1,3,2)
pzplot(cloop2);
title('m. fechada CL2');
axis equal;
subplot(1,3,3)
pzplot(cloop3);
title('m. fechada CL3');
axis equal;
%
% Grafico da resposta a degrau
%
tfinal = 20;
t = 0:0.02:tfinal;
figure(3);
stepplot(cloop1,cloop2,cloop3,t);
grid on
%
% Caracteristicas da resposta a degrau
% Tr - Tempo de subida, Ts - Tempo de acomodacao, Mp - Maximo sobresinal
%
display('caracteristicas da resposta a degrau cloop1')
info=stepinfo(cloop1)
display('caracteristicas da resposta a degrau cloop2')
info=stepinfo(cloop2)
display('caracteristicas da resposta a degrau cloop3')
info=stepinfo(cloop3)
%
% Funcao de transferencia para calculo do esforco de controle u(t)
%
esforco1 =feedback(Hp1,Gp)
esforco2 =feedback(Hp2,Gp)
esforco3 =feedback(Hp3,Gp)
figure(4)
stepplot(esforco1,esforco2,esforco3,t);
grid on
%
% Caracteristicas do esforco de controle
%
display('caracteristicas do esforco de controle de cloop1')
info=stepinfo(esforco1)
display('caracteristicas do esforco de controle de cloop2')
info=stepinfo(esforco2)
display('caracteristicas do esforco de controle de cloop3')
info=stepinfo(esforco3)
