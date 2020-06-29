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

%%% Ganhos de Controle por alocação
%%% Polos desejados
polos = [-4+6*1i, -4-6*1i, -5, -10, -15, -20]; % Primeira Alocacao
polos_v2 = [-2+6*1i, -2-6*1i, -3+4*1i, -3-4*1i, -8, -10]; % Segunda
polos_v3 = [-5.5, -7.5, -9.5, -11.5, -13.5, -15.5]; %terceira
polos_v4 = [-1, -2, -3, -14, -16, -18]; %terceira

%%% Alocação
K = place(A,B,polos);
K_v2 = place(A,B,polos_v2);
K_v3 = place(A,B,polos_v3);
disp(K_v3)
K_v4 = place(A,B,polos_v4);
F = A - B*K;
F_v2 = A - B*K_v2;
F_v3 = A - B*K_v3;
F_v4 = A - B*K_v4;
ee_fechada = ss(F,B,C,D); % Espaco de estados por alocacao
ee_fechada_v2 = ss(F_v2,B,C,D); % Segundo Espaco
ee_fechada_v3 = ss(F_v3,B,C,D); %Terceiro
ee_fechada_v4 = ss(F_v4,B,C,D); %Terceiro

%%% Controle LQ
Q1 = diag([2000, 50, 100, 0.7, 10, 0.7]); %Ref no relatorio
R1 = 0.005;
Q2 = diag([1000, 0.3, 1000, 0.3, 1000, 0.3]); %Ref no relatorio
R2 = 0.2;
Klq1 = lqr(ee,Q1,R1);
Flq1 = A - B*Klq1;
Klq2 = lqr(ee,Q2,R2);
Flq2 = A - B*Klq2;
ee_LQ1 = ss(Flq1,B,C,D); % Espaco de estados por LQ
ee_LQ2 = ss(Flq2,B,C,D); % Espaco de estados por LQ

%%% Matrizes de Transição
Dt = 0.005;
Phi3 = expm(F_v3*Dt);
Phi4 = expm(F_v4*Dt);
Phi_LQ1 = expm(Flq1*Dt);
Phi_LQ2 = expm(Flq2*Dt);
ti = 0;
tf = 1.5;
t = ti:Dt:tf;
x = zeros(size(t));
disp(x)

%%% Estado sem observador
x(1,1) = 0;
x(2,1) = 0;
x(3,1) = 0.25;
x(4,1) = 0;
x(5,1) = 0.25;
x(6,1) = 0;


x3 = x;
x4 = x;
x_LQ1 = x;
x_LQ2 = x;
for i=1:(tf/Dt)
    x3(:,i+1) = Phi3*x3(:,i);
    x_LQ1(:,i+1) = Phi_LQ1*x_LQ1(:,i);
    x4(:,i+1) = Phi4*x4(:,i);
    x_LQ2(:,i+1) = Phi_LQ2*x_LQ2(:,i);
end

%%% Sem observador
figure
plot(t,x3(1,:),t,x3(3,:),t,x3(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
ylim([-0.05 0.3])
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('x (m)','$\theta_1$ (rad)','$\theta_2$ (rad)')

figure
plot(t,x3(2,:),t,x3(4,:),t,x3(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
ylim([-1.2 0.3])
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')

% figure
% plot(t,x4(1,:),t,x4(3,:),t,x4(5,:))
% grid on
% xlabel('Tempo (s)')
% ylabel('Saída no Tempo')
% title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos v4')
% legend('x','$\theta_1$','$\theta_2$')

figure
plot(t,x_LQ1(1,:),t,x_LQ1(3,:),t,x_LQ1(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
ylim([-0.05 0.3])
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('x (m)','$\theta_1$ (rad)','$\theta_2$ (rad)')

figure
plot(t,x_LQ1(2,:),t,x_LQ1(4,:),t,x_LQ1(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
ylim([-1.2 0.3])
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')
 
% figure
% plot(t,x_LQ2(2,:),t,x_LQ2(4,:),t,x_LQ2(6,:))
% grid on
% xlabel('Tempo (s)')
% ylabel('Saída no Tempo')
% title('Simulação da Matriz de Transição de Malha Fechada - Método LQ v2')
% legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')

matriz_latex(Phi_LQ1,'PhiLQ.txt')
matriz_latex(Klq1,'K_LQ.txt')
matriz_latex(Q1,'Q1.txt')

% %DESCOMENTAR PARA VER O VIDEO DO LQ
% x = x_LQ;

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
