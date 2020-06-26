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
%Btil = [0;L;0;0;0;0];

A = Minv*Atil;

B = Minv*Btil;

%C = [1,0,0,0,0,0];
C = eye(6);

D = 0;

ee = ss(A,B,C,D); % Espaço de Estados de malha aberta

%%% Ganhos de Controle por alocação
%%% Polos desejados
polos = [-4+6*1i, -4-6*1i, -5, -10, -15, -20]; % Primeira Alocacao
polos_v2 = [-2+6*1i, -2-6*1i, -3+4*1i, -3-4*1i, -8, -10]; % Segunda

%%% Alocação
K = place(A,B,polos);
K_v2 = place(A,B,polos_v2);
F = A - B*K;
F_v2 = A - B*K_v2;
ee_fechada = ss(F,B,C,D); % Espaco de estados por alocacao
ee_fechada_v2 = ss(F_v2,B,C,D); % Segundo Espaco

% %%% Controle LQ
% Q = diag([15, 15, 15, 15, 15, 15]); %Ref no relatorio
% R = 1;
% Klq = lqr(ee,Q,R);
% Flq = A - B*Klq;
% ee_LQ = ss(Flq,B,C,D); % Espaco de estados por LQ

%%% Alocacao do Observador
polos_ob = [-21, -26, -31, -36, -41, -46]; % Primeira alocacao observ
polos_ob_v2 = [-6+4*1i, -6-4*1i, -8+8*1i, -8-8*1i, -10+2*1i, -10-2*1i]; % Segunda
Ko = (place(A',C',polos_ob))';
Ko_v2 = (place(A',C',polos_ob_v2))';

% %%% Observador LQ
% Qo = diag([5, 5, 5, 5, 5, 5]);
% Ro = 1;
% Ko_lq = (lqr(A',C',Qo,Ro))';

%%% Composicao de Lambda
L11 = A - B*K;
L12 = B*K;
L21 = zeros(6);
L22 = A - Ko*C;
La = [ L11 , L12 ; L21 , L22 ] ;

%%% Lambda Segunda Alocacao
L11_v2 = A - B*K_v2;
L12_v2 = B*K_v2;
L21_v2 = zeros(6);
L22_v2 = A - Ko_v2*C;
La_v2 = [L11_v2 L12_v2; L21_v2 L22_v2];

% %%% Lambda LQ
% L11_lq = A - B*Klq;
% L12_lq = B*Klq;
% L21_lq = L21;
% L22_lq = A - Ko_lq*C;
% L_lq = [L11_lq, L12_lq; L21_lq, L22_lq];


%%% Matrizes de Transição
Dt = 0.01;
Phi = expm(F*Dt);
%Phi_LQ = expm(Flq*Dt);
Phi_lambda = expm(La*Dt);
%Phi_lambda_lq = expm(L_lq*Dt);
Phi_lambda_v2 = expm(La_v2*Dt);
ti = 0;
tf = 20;
t = ti:Dt:tf;
x = zeros(size(t));
xo = x;
%%% Estado com observador
x(1,1) = 0;
x(2,1) = 0;
x(3,1) = 0.1;
x(4,1) = 0;
x(5,1) = 0.2;
x(6,1) = 0;
%%% Estado com observador
xo(1,1) = 0;
xo(2,1) = 0;
xo(3,1) = 0.1;
xo(4,1) = 0;
xo(5,1) = 0.2;
xo(6,1) = 0;
%%% Observador
xo(7,1) = x(1,1) - 0;
xo(8,1) = x(2,1) - 0;
xo(9,1) = x(3,1) - 0;
xo(10,1) = x(4,1) - 0;
xo(11,1) = x(5,1) - 0;
xo(12,1) = x(6,1) - 0;

x_LQ = x;
xo_LQ = xo;
xo_v2 = xo;

for i=1:(tf/Dt)
    x(:,i+1) = Phi*x(:,i);
    %x_LQ(:,i+1) = Phi_LQ*x_LQ(:,i);
    xo(:,i+1) = Phi_lambda*xo(:,i);
    %xo_LQ(:,i+1) = Phi_lambda_lq*xo_LQ(:,i);
    xo_v2(:,i+1) = Phi_lambda_v2*xo_v2(:,i);
end

%%% Plots das matrizes de transição

%%% Com Observ v1
figure
plot(t,xo(1:6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('x','$\dot{x}$','$\theta_1$','$\dot{\theta_1}$','$\theta_2$','$\dot{\theta_2}$')

figure
plot(t,xo(7:12,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Observador - Alocação de Polos')
legend('Erro em x','Erro em $\dot{x}$','Erro em $\theta_1$','Erro em $\dot{\theta_1}$','Erro em $\theta_2$','Erro em $\dot{\theta_2}$')

% %%% Com Observ LQ
% figure
% plot(t,xo_LQ(1:6,:))
% grid on
% xlabel('Tempo (s)')
% ylabel('Saída no Tempo')
% title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
% legend('x','$\dot{x}$','$\theta_1$','$\dot{\theta_1}$','$\theta_2$','$\dot{\theta_2}$')
% 
% figure
% plot(t,xo_LQ(7:12,:))
% grid on
% xlabel('Tempo (s)')
% ylabel('Saída no Tempo')
% title('Simulação da Matriz de Transição de Observador - Método LQ')
% legend('Erro em x','Erro em $\dot{x}$','Erro em $\theta_1$','Erro em $\dot{\theta_1}$','Erro em $\theta_2$','Erro em $\dot{\theta_2}$')

%%% Com Observ v2
figure
plot(t,xo_v2(1:6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos v2')
legend('x','$\dot{x}$','$\theta_1$','$\dot{\theta_1}$','$\theta_2$','$\dot{\theta_2}$')

figure
plot(t,xo_v2(7:12,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Observador - Alocação de Polos v2')
legend('Erro em x','Erro em $\dot{x}$','Erro em $\theta_1$','Erro em $\dot{\theta_1}$','Erro em $\theta_2$','Erro em $\dot{\theta_2}$')

%%% Sem observador
figure
plot(t,x)
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('x','$\dot{x}$','$\theta_1$','$\dot{\theta_1}$','$\theta_2$','$\dot{\theta_2}$')

figure
plot(t,x_LQ)
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Método LQ')
legend('x','$\dot{x}$','$\theta_1$','$\dot{\theta_1}$','$\theta_2$','$\dot{\theta_2}$')

%%% Salva as matrizes de Transição em formato LaTeX
matriz_latex(Phi,'Phi.txt')
% matriz_latex(Phi_LQ, 'Phi_LQ.txt')
% matriz_latex(Phi_lambda, 'Phi_lambda.txt')
% matriz_latex(Phi_lambda_lq, 'Phi_lambda_lq.txt')
% matriz_latex(Phi_lambda_v2, 'Phi_lambda_v2.txt')

% %DESCOMENTAR PARA VER O VIDEO DO LQ
% x = x_LQ;
% %DESCOMENTAR PARA VER O VIDEO DO OBSERVADOR 1
% x = xo;
% %DESCOMENTAR PARA VER O VIDEO DO OBSERVADOR LQ
% x = xo_LQ;
% %DESCOMENTAR PARA VER O VIDEO DO OBSERVADOR 2
% x = xo_v2;

%%% Vídeo da Animação
xA = x(1,:);
yA = zeros(1,length(xA));
yA = yA + 2;
xB = xA + L*sin(x(3,:));
yB = yA - L*cos(x(3,:));
xC = xB + L*sin(x(5,:));
yC = yB - L*cos(x(5,:));
n=length(xA);
v = VideoWriter('C:\Users\Lucas Paiva\Documents\USP\Trabalhos\Controle\T2\Imagens\GIF.avi');
v.Quality = 95;
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
