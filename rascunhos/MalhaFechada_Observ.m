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

%%% Ganhos de Controle por alocação
%%% Polos desejados
polos = [-5.5, -7.5, -9.5, -11.5, -13.5, -15.5];


%%% Alocação
K = place(A,B,polos);
F = A - B*K;
ee_fechada = ss(F,B,C,D); % Espaco de estados por alocacao

%%% Metodo LQ
Q = diag([2000, 50, 100, 0.7, 10, 0.7]);
R = diag([0.005 0.005 0.005]);
Klq = lqr(A,B,Q,R);
Flq = A - B*Klq;
ee_LQ = ss(Flq,B,C,D); % Espaco de estados por LQ

%%% Observador AP
polos_ob = [-16, -18, -20, -22, -24, -26];
Ko = (place(A',C',polos_ob))';

%%% Observador LQ
Qo = diag([5000, 100, 200, 0.07, 50, 0.07]);
Ro = diag([5e-4 5e-4 5e-4]);
Ko_lq = (lqr(A',C',Qo,Ro))';

%%% Composicao de Lambda AP
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
Dt = 0.005;
Phi_lambda = expm(La*Dt);
Phi_lambda_lq = expm(L_lq*Dt);
ti = 0;
tf = 2;
t = ti:Dt:tf;
x = zeros(size(t));

%%% Estado
x(1,1) = 0.2;
x(2,1) = 0;
x(3,1) = 0.25;
x(4,1) = 0;
x(5,1) = 0.25;
x(6,1) = 0;

%%% Observador
x(7,1) = 0;
x(8,1) = x(2,1) - 0;
x(9,1) = x(3,1) - 0;
x(10,1) = 0;
x(11,1) = x(5,1) - 0;
x(12,1) = 0;

x_LQ = x;

for i=1:(tf/Dt)
    x(:,i+1) = Phi_lambda*x(:,i);
    x_LQ(:,i+1) = Phi_lambda_lq*x_LQ(:,i);
end

%%% Plots

%%% Alocação

%%% Estado
figure
plot(t,x(1,:),t,x(3,:),t,x(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('x','$\theta_1$','$\theta_2$')

figure
plot(t,x(2,:),t,x(4,:),t,x(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')

%%% Observador
figure
plot(t,x(7,:),t,x(9,:),t,x(11,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('Erro em x','Erro em $\theta_1$','Erro em $\theta_2$')

figure
plot(t,x(8,:),t,x(10,:),t,x(12,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('Erro em $\dot{x}$','Erro em $\dot{\theta_1}$','Erro em $\dot{\theta_2}$')

%%% LQ

%%% Estado
figure
plot(t,x_LQ(1,:),t,x_LQ(3,:),t,x_LQ(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('x','$\theta_1$','$\theta_2$')

figure
plot(t,x_LQ(2,:),t,x_LQ(4,:),t,x_LQ(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')

%%% Observador
figure
plot(t,x_LQ(7,:),t,x_LQ(9,:),t,x_LQ(11,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('Erro em x','Erro em $\theta_1$','Erro em $\theta_2$')

figure
plot(t,x_LQ(8,:),t,x_LQ(10,:),t,x_LQ(12,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('Erro em $\dot{x}$','Erro em $\dot{\theta_1}$','Erro em $\dot{\theta_2}$')

%%% Descomentar para ver o video do LQ
% x = x_LQ;

%%% Vídeo da Animação
figure
xA = x(1,:);
yA = zeros(1,length(xA));
yA = yA + 2;
xB = xA + L*sin(x(3,:));
yB = yA - L*cos(x(3,:));
xC = xB + L*sin(x(5,:));
yC = yB - L*cos(x(5,:));
n=length(xA);
v = VideoWriter('C:\Users\Lucas Paiva\Documents\USP\Trabalhos\Controle\T2\Imagens\GIF_OBS.avi');
v.Quality = 95;
v.FrameRate = 1/Dt;
open(v);
for i=1:n
Pxa = [ xA(i), xB(i) ];
Pya = [ yA(i), yB(i) ];
Pxb = [ xB(i), xC(i) ];
Pyb = [ yB(i), yC(i) ];
PC1 = [ xA(i)+0.15 , xA(i)+0.15];
PC2 = [ yA(i)+0.03 , yA(i)-0.03];
PC3 = [ xA(i)+0.15 , xA(i)-0.15];
PC4 = [ yA(i)+0.03 , yA(i)+0.03];
PC5 = [ xA(i)-0.15 , xA(i)-0.15];
PC6 = [ yA(i)+0.03 , yA(i)-0.03];
PC7 = [ xA(i)+0.15 , xA(i)-0.15];
PC8 = [ yA(i)-0.03 , yA(i)-0.03];
plot(xA, yA, xA(i),yA(i),'o', xB, yB, xB(i),yB(i),'o', xC, yC, xC(i),yC(i),'o')
axis([-0.6 0.6 0.9 2.1])
grid on
line(Pxa,Pya)
line(Pxb,Pyb)
line(PC1, PC2)
line(PC3, PC4)
line(PC5, PC6)
line(PC7, PC8)
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v)