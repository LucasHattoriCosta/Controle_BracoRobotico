clear; close all; clc

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
FT1 = fts(6,3); % FT relacionando thetadot2 x T2

%%% Consertando FT1

[num,den]=tfdata(FT1,'v');
num2 = [num 0];
den2 = [den 0];
FT1 = tf(num2,den2);

% figure
% margin(FT1)
% legend('Ft original')
% 
% C_PI = tf([63.049 63.049*5.643], [1 0]);
% % C_PID = tf([0.86076 0.86076*57.16 0.86076*908.6], [1 0]);
% C_PID = tf([0.86076/200 0.86076*57.16 0.86076*908.6], [1 0]);
% 
% FTMF = FT1*C_PID;
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
[mag,phase,wout] = bode(FT1);
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

%%% Teste
figure
margin(FT1*Gc)
hold
bode(FT1)
legend('Com Compensador de Bode','Sem Compensador')
 
% %% item d)
% polos = pole(FTMF);
% zeros = zero(FTMF);
% polos_e_zeros = [polos;zeros];
% 
% %controlSystemDesigner(FTMF)
% 
% %%% Polo escolhido iterativamente
% p_det = -0.5e-4+7*1i;
% Im = imag(p_det);
% R = real(p_det);
% 
% figure
% scatter(real(polos_e_zeros), imag(polos_e_zeros))
% grid on
% legend('polos e zeros')
% 
% wn = sqrt(Im^2 + R^2);
% zeta = abs(R/wn);
% 
% thetas = [0 0 0 0 0 0 0 0 0 0 0 0 0];
% imaginarios = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% reais = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% for i=1:14
%     aux = atand((abs(Im - imag(polos_e_zeros(i)))/(abs(R - real(polos_e_zeros(i))))));
%     if real(polos_e_zeros(i)) - R > 0
%         theta = 180 - aux;
%         thetas(i) = theta;
%     else 
%         theta = aux;
%         thetas(i) = theta;
%     end
%     reais(i) = real(polos_e_zeros(i));
%     imaginarios(i) = imag(polos_e_zeros(i));
% end
% reais(6) = real(p_det);
% imaginarios(6) = imag(p_det);
% %thetas(2) = -thetas(2);
% 
% % figure
% % scatter(reais, imaginarios)
% % grid on
% % axis([-6 6 -6 6])
% 
% %%% Determinação do avanço
% phi_LR = - 180 + thetas(1) + thetas(2) + thetas(3) + thetas(4) - thetas(5) - thetas(6) - thetas(7) + thetas(8) + thetas(9) + thetas(10) + thetas(11) + thetas(12) + thetas(13) + thetas(14);
% 
% %%% Determinação do polo e zero do compensador
% bis = (thetas(1))/2;
% 
% xp = abs(Im)*tand(bis - atand(abs(R)/abs(Im)) + phi_LR/2);
% xz = abs(Im)*tand(bis - atand(abs(R)/abs(Im)) - phi_LR/2);
% 
% p_c = xp + abs(R);
% z_c = xz + abs(R);
% 
% %%% Determinação de Kc
% num_aux = [1 z_c];
% den_aux = [1 p_c];
% Gc_aux = tf(num_aux,den_aux);
% G_aux2 = Gc_aux*(FTMF);
% Gc_aux3 = 1/G_aux2;
% K_aux = abs(evalfr(Gc_aux3,p_det));
% Kc = K_aux/0.03134;
% 
% %%% Gc real
% num_Gc_LR = Kc*[1 z_c];
% den_Gc_LR = [1 p_c];
% Gc_LR = tf(num_Gc_LR,den_Gc_LR);
%  
% figure
% margin(Gc_LR*FTMF)
% hold
% bode(FTMF)
% legend('Com Compensador LR','Sem Compensador')
% 
% figure
% rlocus(Gc_LR*FTMF)
% hold
% scatter(R, Im)
% title('Verificação do polo forçado')
% 
% %% item e)
% %%% FTMF LR
% FTMFc_LR = feedback(Gc_LR*FTMF,1);
% FTMFc_bode = feedback(Gc*FTMF,1);
%  
% figure
% step(FTMFc_bode)
% hold
% step(FTMFc_LR)
% legend('comp_bode', 'comp_LR')
% 
% stepinfo(FTMFc_bode)
% stepinfo(FTMFc_LR)
% 
% %% item g)
% figure
% nyquist(FTMF)
% title('original')
% 
% figure
% nyquist(Gc_LR*FTMF)
% title('comp_LR')
% 
% figure
% nyquist(Gc*FTMF)
% title('comp_bode')
% 
% % %% item h)
% % 
% % Kd = 10;
% % Kp = 15;
% % Ki = 35;
% % Gc_PID = tf([Kd Kp Ki],[1 0]);
% % 
% % FT_ITAE = Gc_PID*Gc_LR*G;
% % 
% % FTMF_ITAE = feedback(FT_ITAE,1);
% % 
% % wn_ITAE = (2.163e04)^(1/6);
% % Coef = [1 2.152*wn_ITAE 5.629*wn_ITAE^2 6.934*wn_ITAE^3 6.792*wn_ITAE^4 3.740*wn_ITAE^5  wn_ITAE^6];
% %% item i)
% 
% %controlSystemDesigner(Gc_LR*G)
% 
% %% item j)
% % figure
% % step(FTMF_ITAE)
% 
% %%% step do sistema controlado por LR foi obtido diretamente do
% %%% controlSystemDesigner
% 
% %% item k)
% G_PID_LR = tf([17.142 17.142*2 17.142*2],[1 0]);
% FTMF_PID_LR = feedback(G_PID_LR*Gc*FTMF,1);
% 
% % figure
% % hold on
% % step(FTMF_bode)
% % step(FTMF_LR)
% % step(FTMF_ITAE)
% % step(FTMF_PID_LR)
% % axis([0 1 0 1.8])
% % legend('Bode','LR', 'Controlado com ITAE', 'Controlado com LR')
% 
% %% item l)
% % figure
% % rlocus(FTMF_ITAE)
% % axis([-15 0 -80 80])
% % 
% % figure
% % rlocus(FTMF_ITAE)
% % axis([-1.8 0 -4 4])
% 
% %% item m)
% % nyquist(FTMF_ITAE)
% 
% %% item n)
% %margin(FTMF_ITAE)

nyquist(FT1)