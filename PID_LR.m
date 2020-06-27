%% Start

clear; close all; clc;

%%% Deixa os eixos em LaTeX
set(groot, 'defaultTextInterpreter','latex');

%% Definindo FT's

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

C = eye(6);

D = zeros(6,3);

ee = ss(A,B,C,D); % Espaco de Estados de malha aberta

fts = tf(ee); % Mudanca para FTs

FT_T2_theta2dot = fts(6,3); % FT relacionando thetadot2 x T2

%%% Consertando FT1

[num,den]=tfdata(FT_T2_theta2dot,'v');
num2 = [num 0];
den2 = [den 0];
FT_T2_theta2dot = tf(num2,den2);

%% Ganhos calculados no controlSystemDesigner

C_P = 1.6025;
C_PI = tf([63.049 63.049*5.643], [1 0]);
C_PID = tf([0.86076 0.86076*57.16 0.86076*908.6], [1 0]);

%% FTMFs

FTMF_P = feedback(C_P*FT_T2_theta2dot,1);
FTMF_PI = feedback(C_PI*FT_T2_theta2dot,1);
FTMF_PID = feedback(C_PID*FT_T2_theta2dot,1);

%% Plots
opt = stepDataOptions('StepAmplitude',0.1);
stepplot(FTMF_PI,FTMF_PID,opt)
title('Comparativo das respostas ao degrau - PID por LR')
xlabel('Tempo (s)')
ylabel('Velocidade Angular (rad/s)')
legend('PI','PID')

%% Parametros do Step
stepinfo(FTMF_PI)
stepinfo(FTMF_PID)




