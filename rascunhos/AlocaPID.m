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

%% Ganhos calculados por Alocação de Polos

%%% Controlador PID
Kp_PID = 50.77;
Ki_PID = 136.3;
Kd_PID = 1.31;
C_PID = tf([Kd_PID Kp_PID Ki_PID],[1 0]);
 
%% FTMF
FTMF_PID = feedback(C_PID*FT_T2_theta2dot,1);

%% Plots
figure
opt = stepDataOptions('StepAmplitude',1);
step(FTMF_PID,opt)
title({'Resposta ao degrau da quinta iteração - PID por AP'}, 'Fontsize', 16)
xlabel('Tempo')
ylabel('Velocidade Angular (rad/s)')

%% Parametros do Step
stepinfo(FTMF_PID)
[s,t] = step(FTMF_PID,opt);
RP = 1-s(end);
disp(strcat('RP = ',num2str(RP)))

%% Plots Juntos
Kp_ar = [22.38 15.82 28.13 26.99 50.77 59.29];
Ki_ar = [57.62 47.66 67.92 63.64 136.3 130.86];
Kd_ar = [1.01 0.53 1.44 1.55 1.31 2.33];

figure
hold on
for i=1:6
    K_aux = tf([Kd_ar(i) Kp_ar(i) Ki_ar(i)],[1 0]);
    FTMF_aux = feedback(K_aux*FT_T2_theta2dot,1);
    step(FTMF_aux,opt)
end
[s2,t2] = step(FTMF_aux,opt);
ref = ones(1,size(t2,1));
plot(t2,ref,'Color',[169 169 169]/255)
legend('Conjunto 1','Conjunto 2','Conjunto 3','Conjunto 4','Conjunto 5','Conjunto 6')
title({'Comparativo das respostas ao degrau - PID por AP'}, 'Fontsize', 16)
xlabel('Tempo')
ylabel('Velocidade Angular (rad/s)')