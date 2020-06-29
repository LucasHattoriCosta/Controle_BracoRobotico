%% Start

clear; close all; clc;

%%% Deixa os eixos em LaTeX
set(groot, 'defaultTextInterpreter','latex');

%% Calculando num e den ITAE

syms s KP KD KI;

KPID_num = (KD*s^2+KP*s+KI);
KPID_den = s;

FT_num = 7.281*s^5 + 0.01002*s^4 + 238.1*s^3 + 0.1695*s^2;
FT_den = s^6 + 0.001017*s^5 + 68.78*s^4 + 0.05918*s^3 + 1168*s^2 + 0.8315*s;


Numerador = FT_num * KPID_num;
expand(Numerador);
collect(Numerador);

F_den = FT_den * KPID_den;
expand(F_den);
Denominador = F_den + Numerador;
expand(Denominador);
collect(Denominador);

coef_den_ITAE = coeffs(collect(Denominador));

%% Depois de um trabalho de corno na mão

A = [7.281 0 0;
    0.501/50 0 7.281;
    238.1 7.281 0.501/50;
    0.339/2 0.501/50 238.1;
    0 238.1 0.339/2;
    0 0.339/2 0
    ]

tempo_acomodacao = 1
wn = 4
B = [0;
    -0.001+2.217*wn;
    -343.9/5+6.745*(wn^2);
    -0.2959/5+9.349*(wn^3);
    -1168+11.58*(wn^4);
    -1.663/2+8.68*(wn^5)
    ]



