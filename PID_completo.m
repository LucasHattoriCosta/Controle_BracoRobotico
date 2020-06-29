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
[num,den]=tfdata(FT_T2_theta2dot,'v');
num2 = [num 0];
den2 = [den 0];
FT_T2_theta2dot = tf(num2,den2);

%% Ganhos calculados no controlSystemDesigner para LR

C_P_LR = 1.6025;
C_PI_LR = tf([63.049 63.049*5.643], [1 0]);
C_PID_LR = tf([0.86076/200 0.86076*57.16 0.86076*908.6], [1 0]);

%% FTMFs

FTMF_P_LR = feedback(C_P_LR*FT_T2_theta2dot,1);
FTMF_PI_LR = feedback(C_PI_LR*FT_T2_theta2dot,1);
FTMF_PID_LR = feedback(C_PID_LR*FT_T2_theta2dot,1);

%% Plots
opt = stepDataOptions('StepAmplitude',1);
stepplot(FTMF_PI_LR,FTMF_PID_LR,opt)
title({'Comparativo das respostas ao degrau - PID por LR'}, 'Fontsize', 16)
xlabel('Tempo (s)')
ylabel('Velocidade Angular (rad/s)')
legend({'PI','PID'}, 'Fontsize', 14)

%% Parametros do Step
stepinfo(FTMF_PI_LR)
stepinfo(FTMF_PID_LR)

%% Ganhos calculados por ITAE

cont = 1;
lastcont = 1;
Kd = 0.0:0.01:0.1;%a
Kp = 55:5:65;%b
Ki = 850:20:1250;%c

for a = 1:length(Kd)
    for b = 1:length(Kp)
        for c = 1:length(Ki)
            Ncc = [Kd(a),Kp(b),Ki(c)];
            Dcc = [1, 0];
            Gcc = tf(Ncc,Dcc);
            GccGH = Gcc*FT1;
            T1 = feedback(GccGH,1);
            T = 1;
            delT = 0.01;
            t = 0:delT:T;
            u = ones(1,length(t));
            y = lsim(T1,u,t);
            ITAE = 0;
            for i = 1:length(t)
                inc = t(i)*abs(1-y(i));
                ITAE = ITAE + inc;
            end
            Results(cont,:) = [Kd(a),Kp(b),Ki(c),ITAE];
            cont = cont +1;
            lastcont = lastcont + 1;
            porcent = 100*cont/(length(Kd)*length(Kp)*length(Ki));
            if lastcont > (length(Kd)*length(Kp)*length(Ki))/100
                disp(vpa(porcent,3));
                lastcont = 0;
            end
        end
    end
end
ResultsOrd = sortrows(Results,4);
Kd = ResultsOrd(1,1);
Kp = ResultsOrd(1,2);
Ki = ResultsOrd(1,3);
PID = pid(Kp,Ki,Kd);
FTMF_ITAE = (feedback(PID*FT1,1));
step(FTMF_ITAE)

%% Ganhos calculados por Alocacao de Polos

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
title({'Resposta ao degrau da quinta itera��o - PID por AP'}, 'Fontsize', 16)
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