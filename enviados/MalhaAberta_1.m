%%% Deixa os eixos em LaTeX
set(groot, 'defaultLegendInterpreter','latex');

%%%Espaço de estados
clear; close all; clc;

m = 2;
Mbase = 6;
L = 0.5;
c = 8.5e-5;
b = 7.12e-3;
g = 9.81;

M = [1, 0, 0, 0, 0, 0; 0, (2*m+Mbase)*L, 0, 3*m*(L^2)/2, 0, m*(L^2)/2; 0, 0, 1, 0, 0, 0; 0, 2*m*L, 0, 3*m*(L^2)/2, 0, 2*m*(L^2)/3; 0, 0, 0, 0, 1, 0; 0, m*L/2, 0, m*(L^2)/6, 0, m*(L^2)/3];

I = eye(6);

Minv = I/M;

Atil = [0, 1, 0, 0, 0, 0; 0, -b*L, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, -3*m*L*g/2, -c, -m*L*g/2, 0; 0, 0, 0, 0, 0, 1; 0, 0, 0, c, -m*L*g/2, c];

Btil = [0, 0, 0; L, 0, 0; 0, 0, 0; 0, 1, 0; 0, 0, 0; 0, 0, 1];

A = Minv*Atil;

B = Minv*Btil;

C = [1, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 0, 1];

D = 0;

ee = ss(A,B,C,D);

%%% Matriz de Controlabilidade
Cont = ctrb(A,B);
Posto = rank(Cont);

%%% Matriz de Observabilidade
Ob = obsv(A,C);
Posto2 = rank(Ob);

%%% Matriz de Transição
Dt = 0.01;
Phi = expm(A*Dt);
ti = 0;
tf = 10;
t = ti:Dt:tf;
x = zeros(size(t));
x(1,1) = 0;
x(2,1) = 0;
x(3,1) = 0.1;
x(4,1) = 0;
x(5,1) = 0.2;
x(6,1) = 0;

for i=1:(tf/Dt)
    x(:,i+1) = Phi*x(:,i);
end

p = pole(ee);

% %%% Plots
% figure
% plot(p,'x')
% axis([-0.001 0 -7 7])
% xlabel('Re')
% ylabel('Im')
% title('Polos do Sistema em Malha Aberta')
% 
% figure
% plot(t,x(1,:),t,x(3,:),t,x(5,:))
% xlabel('Tempo (s)')
% ylabel('Saída no Tempo (m ou rad)')
% title('Simulação da Matriz de Transição de Malha Aberta')
% legend('x','$\theta_1$','$\theta_2$')
% 
% figure
% plot(t,x(2,:),t,x(4,:),t,x(6,:))
% xlabel('Tempo (s)')
% ylabel('Saída no Tempo (m/s ou rad/s)')
% title('Simulação da Matriz de Transição de Malha Aberta')
% legend('$\dot{x}$','$\dot{\theta_1}$','$\dot{\theta_2}$')

% matriz_latex(Phi,'C:\Users\Lucas Paiva\Documents\USP\Trabalhos\Controle\T2\Matlabs\Phi.txt')
% matriz_latex(Cont,'C:\Users\Lucas Paiva\Documents\USP\Trabalhos\Controle\T2\Matlabs\Contr.txt')
matriz_latex(Ob,'C:\Users\Lucas Paiva\Documents\USP\Trabalhos\Controle\T2\T2_controlao\enviados\Obsv.txt')

% figure
% rlocus(ee)
% sgrid
% axis([-0.5 0.5 -8 8])
% title('Lugar das raízes para o sistema')
% xlabel('Re')
% ylabel('Im')
