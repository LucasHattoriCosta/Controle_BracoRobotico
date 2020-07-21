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
Btil = [0, 0, 0, 0, 0, 0; L, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0];

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

%%% Composicao de Lambda
L11 = A - B*K;
L12 = B*K;
L21 = zeros(6);
L22 = A - Ko*C;
La = [ L11 , L12 ; L21 , L22 ] ;

%%% Lambda LQ
L11_lq = A - B*Klq;
L12_lq = B*Klq;
L21_lq = L21;
L22_lq = A - Ko_lq*C;
L_lq = [L11_lq, L12_lq; L21_lq, L22_lq];


%%% Matrizes de Transição
Dt = 0.0005;
Phi_lambda = expm(La*Dt);
Phi_lambda_lq = expm(L_lq*Dt);
ti = 0;
tf = 1.5;
t = ti:Dt:tf;
x = zeros(size(t));
xo = x;
%%% Estado com observador
xo(1,1) = 0;
xo(2,1) = 0.1;
xo(3,1) = 0.25;
xo(4,1) = 0;
xo(5,1) = 0.25;
xo(6,1) = 0;
%%% Observador
xo(7,1) = 0;
xo(8,1) = xo(2,1);
xo(9,1) = xo(3,1);
xo(10,1) = 0;
xo(11,1) = xo(5,1);
xo(12,1) = 0;

xo_LQ = xo;

for i=1:(tf/Dt)
    xo(:,i+1) = Phi_lambda*xo(:,i);
    xo_LQ(:,i+1) = Phi_lambda_lq*xo_LQ(:,i);
end

%%% Acelerações

xdd = diff(xo(2,:))./diff(t(1,:));
theta1dd = diff(xo(4,:))./diff(t(1,:));
theta2dd = diff(xo(6,:))./diff(t(1,:));
figure
plot(t(1:end-1),xdd,t(1:end-1),theta1dd,t(1:end-1),theta2dd)
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('$\ddot{x} \; (m/s^2)$','$\ddot{\theta_1}  \; (rad/s^2)$','$\ddot{\theta_2} \; (rad/s^2)$')

%%% Forças

Flin = ((2*m+Mbase)*L*xdd + b*L*xo(2,1:end-1) + (3*m*(L^2)/2)*theta1dd + (m*(L^2)/2)*theta2dd)/L;
Torq1 = (3*m*(L^2)/2)*theta1dd + (m*(L^2))*theta2dd + 2*m*L*xdd + (3*m*L*g/2)*xo(3,1:end-1) + c*xo(4,1:end-1) + (m*L*g/2)*xo(5,1:end-1);
Torq2 = (m*(L^2)/3)*theta2dd + (m*(L^2)/6)*theta1dd + (m*L/2)*xdd - c*xo(4,1:end-1) + (m*L*g/2)*xo(5,1:end-1) + c*xo(6,1:end-1);

figure
plot(t(1:end-1),Flin,t(1:end-1),Torq1,t(1:end-1),Torq2)
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Força dos Atuadores - Alocação de Polos')
legend('F (N)','$T_1 \; (N \cdot m)$','$T_2 \; (N \cdot m)$')
axis([0 0.5 -80 80])

%%% Potencias

Pot1 = Flin.*xo(2,1:end-1);
Pot2 = Torq1.*xo(4,1:end-1);
Pot3 = Torq2.*xo(6,1:end-1);

figure
plot(t(1:end-1),Pot1,t(1:end-1),Pot2,t(1:end-1),Pot3)
grid on
xlabel('Tempo (s)')
ylabel('Potência (W)')
title('Potência dos Atuadores - Alocação de Polos')
legend('$Pot_F$','$Pot_{T_1}$','$Pot_{T_2}$')
axis([0 0.5 -20 15])
