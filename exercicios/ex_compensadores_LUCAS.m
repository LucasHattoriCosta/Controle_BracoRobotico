clc;close all;

%%% Deixa os eixos em LaTeX
set(groot, 'defaultLegendInterpreter','latex');
%%  Definicao da FT
K = 1;

num = [K K];
den = [1 3 12 -16 0];
G = tf(num, den);
%% Item b
K = 29.5;

% figure
% margin(K * G);
% disp('Polos: ');
% disp(pole(K * G));
% disp(allmargin(K * G));

%% Item c - Para descobrir phi_m

% margem_fase = 0;
% 
% for avanco_alvo = linspace(30,70,50);
%     %%% Passo 1
%     margem_fase_original = 3.1898;
%     %%% Passo 2
%     phi_m = avanco_alvo - margem_fase_original;
%     %%% Passo 3
%     alpha = (sind(phi_m) + 1) / (1 - sind(phi_m));
%     %%% Passo 4
%     g_m = 10 * log10(1/alpha);
%     [mag,phase,wout] = bode(K * G);
%     mag = squeeze(mag);
%     phase = squeeze(phase);
%     w_m = interp1(20*log10(mag), wout, g_m);
%     %%% Passo 5
%     p = sqrt(alpha) * w_m;
%     z = w_m / sqrt(alpha);
%     %%% Passo 6
%     num_Gc = [alpha alpha*z];
%     den_Gc = [1 p];
%     Gc = tf(num_Gc, den_Gc);
%     
%     [gc, phase_c, wcg, wcp] = margin(K * G *Gc);
%     
%     if phase_c > margem_fase
%         margem_fase = phase_c;
%         disp(phi_m);
%         disp(margem_fase);
%     end 
% end

%% Item c - recalcula com a maior margem de fase

% Phi_m determinado
phi_m = 66.8102;
%%% Passo 3
alpha = (sind(phi_m) + 1) / (1 - sind(phi_m));
%%% Passo 4
g_m = 10 * log10(1/alpha);
[mag,phase,wout] = bode(K * G);
mag = squeeze(mag);
phase = squeeze(phase);
w_m = interp1(20*log10(mag), wout, g_m);
%%% Passo 5
p = sqrt(alpha) * w_m;
z = w_m / sqrt(alpha);
%%% Passo 6
Gc_bode = tf([alpha alpha*z], [1 p]);

figure
% margin(K * G * Gc_bode)
% hold
% bode(K * G)
% legend('Com compensador Bode', "Sem compensador")

% %%% FTMF Bode
FTMF_bode = feedback(Gc_bode*K*G,1);

% figure
% step(FTMF_bode)
% hold
% stepinfo(FTMF_bode)
% title('Resposta ao degrau para sistema com compensador - Bode')

%% Item d

% real_pesc = input('Real_pesc :');
% imag_pesc = input('Imag_pesc :');

real_pesc = -0.65;
imag_pesc = 4.5;

pesc = complex(real_pesc, imag_pesc);
p1 = 1;
p2 = 0;
p3 = -1;
p4 = -2 + 3.464*1i;
p5 = -2 - 3.464*1i;
p = [p1 p2 p3 p4 p5];

thetas = [];
for i=1:5
    thetaproto = atand((abs(imag(pesc)-imag(p(i)))/(abs(real(pesc)-real(p(i))))));
    if real(p(i)) - real(pesc) > 0
        theta = 180 - thetaproto;
        thetas = [thetas theta];
    else 
        theta = thetaproto;
        thetas = [thetas theta];
    end
end
aux = (- thetas(1) - thetas(2) + thetas(3) - thetas(4) - thetas(5));
avanco = -180 - aux;
phi_2 = avanco/2;
bissec = thetas(2)/2;

AngZ = bissec - phi_2;
AngP = bissec + phi_2;
Lado = abs(pesc);
Ang0 = 180 - thetas(2);

%%% Da Prime, P e Z vão aparecer em b
Z = Lado * sind(AngZ) / sind(180 - Ang0 - AngZ);
P = Lado * sind(AngP) / sind(180 - Ang0 - AngP);

%%% Achando Kc
Gc_aux = tf([1 Z],[1 P]);
Kc_LR = abs(evalfr(1/(Gc_aux*G),pesc))/K;

%%% Gc LR
Gc_LR = tf([Kc_LR Kc_LR*Z],[1 P]);

% figure
% margin(Gc_LR*K*G)
% hold on
% bode(K*G)
% hold on
% bode(Gc_bode * K * G)
% legend('Com Compensador LR','Sem compensador','Com compensador Bode')

%%% FTMF LR
FTMF_LR = feedback(Gc_LR*K*G,1);

% figure
% step(FTMF_LR)
% stepinfo(FTMF_LR)
% title('Resposta ao degrau para sistema com compensador - LR')

%% Item 3 - steps dos dois compensadores

figure
step(FTMF_bode)
hold on
step(FTMF_LR)
legend({'Compensador via Bode', 'Compensador via LR'},'FontSize',14)
title({'Comparação das respostas ao degrau dos sistemas compensados'}, 'Fontsize',16)

%% Item H - ITAE

syms s KP KD KI;

G_num = (s + 1);
G_den = (s^4 + 3*s^3 + 12*s^2 - 16*s);
KPID_num = (KD*s^2+KP*s+KI);
KPID_den = s;
Gc_num = 6.874*s + 10.23;
Gc_den = s + 13.89;

Numerador_ITAE = K * G_num * Gc_num * KPID_num;
expand(Numerador_ITAE);
collect(Numerador_ITAE);

F_den = G_den * Gc_den * KPID_den;
expand(F_den);
Denominador_ITAE = F_den + Numerador_ITAE;
expand(Denominador_ITAE);
coef_den_ITAE = coeffs(collect(Denominador_ITAE));

wn_ITAE = coef_den_ITAE(13)/2.152;
KD = (5.629 * wn_ITAE^2 - coef_den_ITAE(11)) / coef_den_ITAE(12);
KP = (6.934 * wn_ITAE^3 - coef_den_ITAE(8) - coef_den_ITAE(10) * KD) / coef_den_ITAE(9);
KI = (6.792 * wn_ITAE^4 - coef_den_ITAE(4) - coef_den_ITAE(5) * KP - coef_den_ITAE(7) * KD) / coef_den_ITAE(6);

KD = eval(KD);
KP = eval(KP);
KI = eval(KI);

G_ITAE = tf([KD KP KI],[1 0]);

FTMF_PID_ITAE = feedback(G_ITAE * Gc_LR * K * G, 1);


%% Item I - Criando a FT do compensador via LR e PID via LR

% FTMF_PID_LR = feedback(C_PID_LR * G_PID_LR, 1);
% step(FTMF_PID_LR)

%% Item J - Respostas ao degrau para ambos PID

figure
step(FTMF_PID_ITAE)
hold on
step(FTMF_PID_LR)
legend({'PID via ITAE', 'PID via LR'},'FontSize',14)
title({'Comparação das respostas ao degrau dos sistemas compensados com PID'}, 'Fontsize',16)

%% Item K - Todos os step
close all;

figure
step(FTMF_bode)
hold on
step(FTMF_LR)
hold on
step(FTMF_PID_ITAE)
hold on
step(FTMF_PID_LR)
legend({'Compensador Bode', 'Compensador LR', 'PID ITAE', 'PID LR'},'FontSize',14)
title({'Comparação das respostas ao degrau dos sistemas projetados'}, 'Fontsize',16)

%% Item L - RLocus

rlocus(C_PID_LR * G_PID_LR)
title({'Lugar das raízes do sistema com compensador LR e PID LR'}, 'Fontsize',16)