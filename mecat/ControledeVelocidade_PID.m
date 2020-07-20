%
% Controle de velocidade angular do motor eletrico CC
% Testes com controle PID
%
% Definicao dos valores dos parametros do sistema
%
% Constante do tacometro
Ktac = 0.48;
% Parametros da funcao de transferencia Gw(s)
K1 = 100;
Km = 2.083;
Kg = 0.1;
a  = 100;
am = 1.71;
%
% Funcao de transferencia da velocidade angular do sistema
s = tf('s');
Gw = (K1*Km*Kg)/((s+a)*(s+am))
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
% Controlador 1
Kp1 = 3
Ti1 = Inf;
Td1 = 0;
N1  = Inf;
H1  = pidstd(Kp1,Ti1,Td1,N1)
Hw1 = Ktac*H1
% Controlador 2
Kp2 = 5
Ti2 = Inf;
Td2 = 0;
N2  = Inf;
H2  = pidstd(Kp2,Ti2,Td2,N2)
Hw2 = Ktac*H2
% Controlador 3
Kp3 = 7
Ti3 = Inf;
Td3 = 0;
N3  = Inf;
H3  = pidstd(Kp3,Ti3,Td3,N3)
Hw3 = Ktac*H3
%
% Definicao da malha aberta
% Sao definidos 3 sistemas distintos
%
GHw1 = Hw1*Gw
GHw2 = Hw2*Gw
GHw3 = Hw3*Gw
%
% Polos e zeros de malha aberta
%
display('Polos e zeros - GHw1')
pole(GHw1)
zero(GHw1)
display('Polos e zeros - GHw2')
pole(GHw2)
zero(GHw2)
display('Polos e zeros - GHw3')
pole(GHw3)
zero(GHw3)
%
% Plot com os polos e zeros de malha aberta
%
figure(1);
subplot(1,3,1);
pzplot(GHw1);
title('m. aberta GHw1');
axis equal;
subplot(1,3,2);
pzplot(GHw2);
title('m. aberta GHw2');
axis equal;
subplot(1,3,3);
pzplot(GHw3);
title('m. aberta GHw3');
axis equal;
%
% definicao da malha fechada
% realimentacao unitaria
%
% Funcao de transferencia em malha fechada
% pode ser definida em Matlab utilizando-se
% O comando feedback(S1,S2) 
% Onde o sistema S1 se encontra na malha direta
% e S2 se encontra na malha de realimentacao
% No caso desse sistema a malha direta
% e' G(s)H(s) e malha de realimentacao e' 
% unitaria
%
% R(s)  E(s)|------|  |------|  U(s)
%---->(+)---| H(s) |--| G(s) |------->
%    _ ^    |------|  |------|    |
%      |---------------------------                       
%      
cloop1 = feedback(GHw1,1)
cloop2 = feedback(GHw2,1)
cloop3 = feedback(GHw3,1)
%
% Polos e zeros de malha fechada
% comando pole() obtem os polos do sistema
% e zero() obtem os zeros do sistema
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
% comando pzplot()
% plot os polos e zeros
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
tfinal = 12;
t = 0:0.02:tfinal;
figure(3);
%
% comando stepplot() simula e plota
% a resposta a degrau unitaria
%
stepplot(cloop1,cloop2,cloop3,t);
grid on
%
% Caracteristicas da resposta a degrau
% Tr - Tempo de subida, Ts - Tempo de acomodacao, Mp - Maximo sobresinal
%
% Comando stepinfo()
% calcula as caracteristicas da resposta a degrau
display('caracteristicas da resposta a degrau cloop1')
info=stepinfo(cloop1)
display('caracteristicas da resposta a degrau cloop2')
info=stepinfo(cloop2)
display('caracteristicas da resposta a degrau cloop3')
info=stepinfo(cloop3)
%
% Funcao de transferencia para calculo do esforco de controle u(t)
% o sinal u(t) pode ser calculado definindo-se
% um sistema de controle em malha fechada onde H(s)
% esta na malha direta e G(s) na malha de realimentacao
%
% R(s)  E(s)|------|        U(s)
%---->(+)---| H(s) |------------>
%    _ ^    |------|    |
%      |                |
%      |    |------|    |
%      |----| G(s) |<----
%           |------|
%      
esforco1 = feedback(Hw1,Gw)
esforco2 = feedback(Hw2,Gw)
esforco3 = feedback(Hw3,Gw)
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
