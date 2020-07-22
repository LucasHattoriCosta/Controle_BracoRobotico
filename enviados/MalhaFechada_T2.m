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

%%% Matrizes de Transição
Dt = 0.005;
Phi = expm(F*Dt);
ti = 0;
tf = 1.5;
t = ti:Dt:tf;
x = zeros(size(t));

%%% Estado sem observador
x(1,1) = 0.2;
x(2,1) = 0;
x(3,1) = 0.25;
x(4,1) = 0;
x(5,1) = 0.25;
x(6,1) = 0;

for i=1:(tf/Dt)
    x(:,i+1) = Phi*x(:,i);
end

%%% Acelerações

xdd = diff(x(2,:))./diff(t(1,:));
theta1dd = diff(x(4,:))./diff(t(1,:));
theta2dd = diff(x(6,:))./diff(t(1,:));
figure
plot(t(1:end-1),xdd,t(1:end-1),theta1dd,t(1:end-1),theta2dd)
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo (m/s^2 e rad/s^2)')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('$\ddot{x} \; (m/s^2)$','$\ddot{\theta_1}  \; (rad/s^2)$','$\ddot{\theta_2} \; (rad/s^2)$')

%%% Forças

Flin = ((2*m+Mbase)*L*xdd + b*L*x(2,1:end-1) + (3*m*(L^2)/2)*theta1dd + (m*(L^2)/2)*theta2dd)/L;
Torq1 = (3*m*(L^2)/2)*theta1dd + (m*(L^2))*theta2dd + 2*m*L*xdd + (3*m*L*g/2)*x(3,1:end-1) + c*x(4,1:end-1) + (m*L*g/2)*x(5,1:end-1);
Torq2 = (m*(L^2)/3)*theta2dd + (m*(L^2)/6)*theta1dd + (m*L/2)*xdd - c*x(4,1:end-1) + (m*L*g/2)*x(5,1:end-1) + c*x(6,1:end-1);

figure
plot(t(1:end-1),Flin,t(1:end-1),Torq1,t(1:end-1),Torq2)
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Força dos Atuadores - Alocação de Polos')
legend('F (N)','$T_1 \; (N \cdot m)$','$T_2 \; (N \cdot m)$')

%%% Potencias

Pot1 = Flin.*x(2,1:end-1);
Pot2 = Torq1.*x(4,1:end-1);
Pot3 = Torq2.*x(6,1:end-1);

figure
plot(t(1:end-1),Pot1,t(1:end-1),Pot2,t(1:end-1),Pot3)
grid on
xlabel('Tempo (s)')
ylabel('Potência (W)')
title('Potência dos Atuadores - Alocação de Polos')
legend('$Pot_F$','$Pot_{T_1}$','$Pot_{T_2}$')

%%% Plots

figure
plot(t,x(1,:),t,x(3,:),t,x(5,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('x (m)','$\theta_1$ (rad)','$\theta_2$ (rad)')

figure
plot(t,x(2,:),t,x(4,:),t,x(6,:))
grid on
xlabel('Tempo (s)')
ylabel('Saída no Tempo')
title('Simulação da Matriz de Transição de Malha Fechada - Alocação de Polos')
legend('$\dot{x}$ (m/s)','$\dot{\theta_1}$ (rad/s)','$\dot{\theta_2}$ (rad/s)')
axis([0 1.5 -1.2 0.1])


% %%% Vídeo da Animação
% xA = x(1,:);
% yA = zeros(1,length(xA));
% yA = yA + 2;
% xB = xA + L*sin(x(3,:));
% yB = yA - L*cos(x(3,:));
% xC = xB + L*sin(x(5,:));
% yC = yB - L*cos(x(5,:));
% n=length(xA);
% v = VideoWriter('C:\Users\Lucas Paiva\Documents\USP\Trabalhos\Controle\T2\Imagens\GIF.avi');
% v.Quality = 95;
% v.FrameRate = 1/Dt;
% open(v);
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
