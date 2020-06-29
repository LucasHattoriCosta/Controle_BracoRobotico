clear; close all; clc

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
C = eye(6);
D = zeros(6,3);

ee = ss(A,B,C,D); % Espaço de Estados de malha aberta
fts = tf(ee); % Mudança para FTs
FT1 = fts(6,3); % FT relacionando thetadot2 x T2

%%% Consertando FT1

[num,den]=tfdata(FT1,'v');
num2 = [num 0];
den2 = [den 0];
FT1 = tf(num2,den2);

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