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

C = eye(6);

D = zeros(6,3);

ee = ss(A,B,C,D); % Espaço de Estados de malha aberta

fts = tf(ee); % Mudança para FTs

FT1 = fts(1,1); % FT relacionando F x x

FT2 = fts(3,2); % FT relacionando T1 x theta1

FT3 = fts(5,3); % FT relacionando T2 x theta2

%%% Para FT1, inicialmente


%%% Dimensionando Kp
[r,k] = rlocus(FT1);
kp = [];
for i=1:size(k,2)
    cont = 0;
    for j=1:size(r,1)
        if r(j,i) <= 0
            cont = cont+1;
        end
    end
    if cont == size(r,1)
        kp = [kp, k(i)];
    end
end

%%% Kp escolhido

Kp = kp(end);

%%% Dimensionando Ki

[N, D] = tfdata(FT1,'v');
D_aux = D + Kp*N;
D_2 = [D_aux,0];
FT1_aux = tf(N,D_2);

[r_2,k_2] = rlocus(FT1_aux);
ki = [];
for i=1:size(k_2,2)
    cont = 0;
    for j=1:size(r_2,1)
        if r_2(j,i) <= 0
            cont = cont+1;
        end
    end
    if cont == size(r_2,1)
        ki = [ki, k_2(i)];
    end
end

%%% Ki escolhido
Ki = ki(end);

%%% Dimensionando Kd
N_3 = [N,0,0];
D_3 = conv([D,0] + ([0,0,0,0,0,0,Kp,0] + [0,0,0,0,0,0,0,Ki]),N);
FT2_aux = tf(N_3,D_3);

[r_3,k_3] = rlocus(FT2_aux);
kd = [];
for i=1:size(k_3,2)
    cont = 0;
    for j=1:size(r_3,1)
        if r_3(j,i) <= 0
            cont = cont+1;
        end
    end
    if cont == size(r_3,1)
        kd = [kd, k_3(i)];
    end
end
disp('Kp = ')
disp(kp)
disp('Ki = ')
disp(ki)
disp('Kd = ')
disp(kd)