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
fts = tf(ee); % Mudança para FTs

FT1 = fts(1,1); % FT relacionando F x x
FT2 = fts(3,2); % FT relacionando T1 x theta1
FT3 = fts(5,3); % FT relacionando T2 x theta2

Kd = 2.7413e6;
Kp = 2*Kd;
Ki = 2*Kd;

Num_K = [Kd Kp Ki];
Den_K = [1 0];
K = tf(Num_K,Den_K);

controlSystemDesigner(FT1)

%%% Matrizes de Transição
Dt = 0.0005;
    ti = 0;
tf = 1.5;
t = ti:Dt:tf;
x = zeros(size(t));

% xo_LQ = xo;
% 
% for i=1:(tf/Dt)
%     xo(:,i+1) = Phi_lambda*xo(:,i);
%     xo_LQ(:,i+1) = Phi_lambda_lq*xo_LQ(:,i);
% end


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
