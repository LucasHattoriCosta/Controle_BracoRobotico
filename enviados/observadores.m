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
C = eye(6);
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

%%% Plots das matrizes de transição
%%% Com Observ v1
figure
plot(t,xo(1,:), t,xo(3,:), t,xo(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('x','$\theta_1$','$\theta_2$')
ylim([-0.05 0.3])

figure
plot(t,xo(2,:), t,xo(4,:), t,xo(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')
ylim([-1.2 0.3])

figure
plot(t,xo(8,:), t,xo(9,:), t,xo(11,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Observador - Alocação de Polos')
legend('Erro em $\dot{x}$','Erro em $\theta_1$','Erro em $\theta_2$')
xlim([0 0.05])
ylim([0 0.25])

%%% Com Observ LQ
figure
plot(t,xo_LQ(1,:),t,xo_LQ(3,:),t,xo_LQ(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('x','$\theta_1$','$\theta_2$')

figure
plot(t,xo_LQ(2,:),t,xo_LQ(4,:),t,xo_LQ(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')
ylim([-1.2 0.3])

figure
plot(t,xo_LQ(8,:),t,xo_LQ(9,:),t,xo_LQ(11,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Observador - Método LQ')
legend('Erro em $\dot{x}$','Erro em $\theta_1$','Erro em $\theta_2$')
xlim([0 0.05])
ylim([0 0.25])

%%% Salva as matrizes de Transição em formato LaTeX
% matriz_latex(Phi_lambda, 'Phi_lambda.txt')
% matriz_latex(Phi_lambda_lq, 'Phi_lambda_lq.txt')
% matriz_latex(Phi_lambda_v2, 'Phi_lambda_v2.txt')

% %DESCOMENTAR PARA VER O VIDEO DO OBSERVADOR 1
% x = xo;
% %DESCOMENTAR PARA VER O VIDEO DO OBSERVADOR LQ
% x = xo_LQ;
% %DESCOMENTAR PARA VER O VIDEO DO OBSERVADOR 2
% x = xo_v2;

% %%% Vídeo da Animação
% xA = x(1,:);
% yA = zeros(1,length(xA));
% yA = yA + 2;
% xB = xA + L*sin(x(3,:));
% yB = yA - L*cos(x(3,:));
% xC = xB + L*sin(x(5,:));
% yC = yB - L*cos(x(5,:));
% n=length(xA);
% v = VideoWriter('GIF.avi');
% v.Quality = 95;
% v.FrameRate = 1/Dt;
% open(v);
% figure
% for i=1:n
% Pxa = [ xA(i), xB(i) ];
% Pya = [ yA(i), yB(i) ];
% Pxb = [ xB(i), xC(i) ];
% Pyb = [ yB(i), yC(i) ];
% PC1 = [ xA(i)+0.15 , xA(i)+0.15];
% PC2 = [ yA(i)+0.03 , yA(i)-0.03];
% PC3 = [ xA(i)+0.15 , xA(i)-0.15];
% PC4 = [ yA(i)+0.03 , yA(i)+0.03];
% PC5 = [ xA(i)-0.15 , xA(i)-0.15];
% PC6 = [ yA(i)+0.03 , yA(i)-0.03];
% PC7 = [ xA(i)+0.15 , xA(i)-0.15];
% PC8 = [ yA(i)-0.03 , yA(i)-0.03];
% plot(xA, yA, xA(i),yA(i),'o', xB, yB, xB(i),yB(i),'o', xC, yC, xC(i),yC(i),'o')
% axis([-0.6 0.6 0.9 2.1])
% grid on
% line(Pxa,Pya)
% line(Pxb,Pyb)
% line(PC1, PC2)
% line(PC3, PC4)
% line(PC5, PC6)
% line(PC7, PC8)
% frame = getframe(gcf);
% writeVideo(v,frame);
% end
% close(v)
