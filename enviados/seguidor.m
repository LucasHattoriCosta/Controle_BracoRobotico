clear; close all; clc;

%%% Deixa os eixos em LaTeX
set(groot, 'defaultLegendInterpreter','latex');

%%% Parametros
m = 2;
Mbase = 6;
L = 0.5;
c = 8.5e-5;
b = 7.12e-3;
g = 9.81;

%%% Espaco de estados
M = [1, 0, 0, 0, 0, 0; 0, (2*m+Mbase)*L, 0, 3*m*(L^2)/2, 0, m*(L^2)/2; 0, 0, 1, 0, 0, 0; 0, 2*m*L, 0, 3*m*(L^2)/2, 0, 2*m*(L^2)/3; 0, 0, 0, 0, 1, 0; 0, m*L/2, 0, m*(L^2)/6, 0, m*(L^2)/3];
I = eye(6);
Minv = I/M;
Atil = [0, 1, 0, 0, 0, 0; 0, -b*L, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, -3*m*L*g/2, -c, -m*L*g/2, 0; 0, 0, 0, 0, 0, 1; 0, 0, 0, c, -m*L*g/2, c];
Btil = [0, 0, 0; L, 0, 0; 0, 0, 0; 0, 1, 0; 0, 0, 0; 0, 0, 1];

A = Minv*Atil;
B = Minv*Btil;
C = [1, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 0, 1];
D = 0;

ee = ss(A,B,C,D); % Espaço de Estados de malha aberta

%%% Alocação
%%% Polos desejados
polos = [-5.5, -7.5, -9.5, -11.5, -13.5, -15.5]; %terceira
K = place(A,B,polos);
F = A - B*K;
ee_fechada = ss(F,B,C,D); %Terceiro

%%% Controle LQ
Q = diag([2000, 50, 100, 0.7, 10, 0.7]); %Ref no relatorio
R = 0.005;
Klq = lqr(ee,Q,R);
Flq = A - B*Klq;
ee_LQ = ss(Flq,B,C,D); % Espaco de estados por LQ

%%% Alocacao do Observador
polos_ob = [-55, -75, -95, -115, -135, -155]; % Primeira alocacao observ
Ko = (place(A',C',polos_ob))';

%%% Observador LQ
Qo = diag([1, 100, 100, 1, 100, 1]);
Ro = 0.005;
Ko_lq = (lqr(A',C',Qo,Ro))';

t = 0:0.01:10;
u = [sin(t); sin(t); sin(t)];%zeros(2,size(t,2))];

R_ap = reg(ee,K,Ko);
FTs_ap = tf(R_ap);
R_lq = reg(ee,Klq,Ko_lq);
FTs_lq = tf(R_lq);

% figure
% lsim(-R_lq(3,3),u,t)
% title('Resposta do sistema ao seguidor senoidal - método LQ')
% xlabel('Tempo')
% ylabel('Velocidade angular (rad/s)')

figure
step(-R_lq(3,3))
title('Resposta do sistema ao seguidor em degrau unitário - método LQ')
xlabel('Tempo')
ylabel('Velocidade angular (rad/s)')