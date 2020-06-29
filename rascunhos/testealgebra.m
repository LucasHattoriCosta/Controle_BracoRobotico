syms s KP KD KI;


K = 29.5;
G_num = (s + 1);
G_den = (s^4 + 3*s^3 + 12*s^2 - 16*s);
KPID_num = (KD*s^2+KP*s+KI);
KPID_den = s;
Gc_num = 6.875*s + 10.23;
Gc_den = s + 13.89; 
Numerador = K * G_num * Gc_num * KPID_num;
expand(Numerador);
collect(Numerador)
F_den = G_den * Gc_den * KPID_den;
expand(F_den);
Denominador =F_den + Numerador;
expand(Denominador);
collect(Denominador)
