%%% Exercício Compensadores
clc; clear; close all;

%% item b)
K = 29.5;
%K = 17.5;
%K = 1;
numerador = K*[1 1];
denominador = [1 3 12 -16 0];

G = tf(numerador, denominador);
G_semK = tf([1 1],denominador);

%figure
%margin(G)

%% item c)

%%% Passo 1
Pm = 3.19;
%%% Passo 2
phi_m = 60 + 5;
%phi_m = 40 + 5;
%%% Passo 3
alpha = (sind(phi_m) + 1)/(1 - sind(phi_m));
%%% Passo 4
g_m = 10 * log10(1/alpha);
[mag,phase,wout] = bode(G);
mag = squeeze(mag);
phase = squeeze(phase);
w_m = interp1(20*log10(mag), wout, g_m);
%%% Passo 5
z = sqrt((w_m^2)/alpha);
p = z*alpha;

%%% Compensador
num_Gc = alpha*[1 z];
den_Gc = [1 p];
Gc = tf(num_Gc,den_Gc);

% %%% Teste
% figure
% margin(G*Gc)
% hold
% bode(G)
% legend('Com Compensador','Sem Compensador')
 
%% item d)
polos = pole(G);
polos_e_zero = [polos;-1];

%%% Polo escolhido iterativamente
p_det = -0.5+5.3*1i;
Im = imag(p_det);
R = real(p_det);

wn = sqrt(Im^2 + R^2);
zeta = abs(R/wn);

thetas = [0 0 0 0 0];
imaginarios = [0 0 0 0 0 0];
reais = [0 0 0 0 0 0];
for i=1:5
    aux = atand((abs(Im - imag(polos_e_zero(i)))/(abs(R - real(polos_e_zero(i))))));
    if real(polos_e_zero(i)) - R > 0
        theta = 180 - aux;
        thetas(i) = theta;
    else 
        theta = aux;
        thetas(i) = theta;
    end
    reais(i) = real(polos_e_zero(i));
    imaginarios(i) = imag(polos_e_zero(i));
end
reais(6) = real(p_det);
imaginarios(6) = imag(p_det);
%thetas(2) = -thetas(2);

% figure
% scatter(reais, imaginarios)
% grid on
% axis([-6 6 -6 6])

%%% Determinação do avanço
phi_LR = - 180 + thetas(1) + thetas(2) + thetas(3) + thetas(4) - thetas(5);

%%% Determinação do polo e zero do compensador
bis = (thetas(1))/2;

xp = abs(Im)*tand(bis - atand(abs(R)/abs(Im)) + phi_LR/2);
xz = abs(Im)*tand(bis - atand(abs(R)/abs(Im)) - phi_LR/2);

p_c = xp + abs(R);
z_c = xz + abs(R);

%%% Determinação de Kc
num_aux = [1 z_c];
den_aux = [1 p_c];
Gc_aux = tf(num_aux,den_aux);
G_aux2 = Gc_aux*(G_semK);
Gc_aux3 = 1/G_aux2;
K_aux = abs(evalfr(Gc_aux3,p_det));
Kc = K_aux/K;

%%% Gc real
num_Gc_LR = Kc*[1 z_c];
den_Gc_LR = [1 p_c];
Gc_LR = tf(num_Gc_LR,den_Gc_LR);
 
% figure
% margin(Gc_bode*G)
% hold
% bode(G)
% legend('Com Compensador','Sem Compensador')

% figure
% rlocus(Gc_bode*G)
% hold
% scatter(R, Im)
% axis([-3 2 -6 6])

%% item e)
%%% FTMF LR
FTMF_LR = feedback(Gc_LR*G,1);
FTMF_bode = feedback(Gc*G,1);
 
% figure
% step(FTMF_bode)
% stepinfo(FTMF_bode)

% figure
% step(FTMF_LR)
% stepinfo(FTMF_LR)

%% item g)
% figure
% nyquist(G)
% axis([-8 0 -4 4])

% figure
% nyqlog(Gc_LR*G)
% axis([-10 0 -5 5])
% 
% figure
% nyqlog(Gc*G)
% axis([-10 0 -5 5])

%% item h)

Kd = 10;
Kp = 15;
Ki = 35;
Gc_PID = tf([Kd Kp Ki],[1 0]);

FT_ITAE = Gc_PID*Gc_LR*G;

FTMF_ITAE = feedback(FT_ITAE,1);

wn_ITAE = (2.163e04)^(1/6);
Coef = [1 2.152*wn_ITAE 5.629*wn_ITAE^2 6.934*wn_ITAE^3 6.792*wn_ITAE^4 3.740*wn_ITAE^5  wn_ITAE^6];
%% item i)

%controlSystemDesigner(Gc_LR*G)

%% item j)
% figure
% step(FTMF_ITAE)

%%% step do sistema controlado por LR foi obtido diretamente do
%%% controlSystemDesigner

%% item k)
G_PID_LR = tf([17.142 17.142*2 17.142*2],[1 0]);
FTMF_PID_LR = feedback(G_PID_LR*Gc*G,1);

% figure
% hold on
% step(FTMF_bode)
% step(FTMF_LR)
% step(FTMF_ITAE)
% step(FTMF_PID_LR)
% axis([0 1 0 1.8])
% legend('Bode','LR', 'Controlado com ITAE', 'Controlado com LR')

%% item l)
% figure
% rlocus(FTMF_ITAE)
% axis([-15 0 -80 80])
% 
% figure
% rlocus(FTMF_ITAE)
% axis([-1.8 0 -4 4])

%% item m)
nyquist(FTMF_ITAE)

%% item n)
%margin(FTMF_ITAE)