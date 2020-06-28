import numpy as np

A = np.array([[7.281, 0, 0],
 [0.501/50, 0, 7.281],
 [238.1, 7.281, 0.501/50],
 [0.339/2, 0.501/50, 238.1],
 [0,238.1,0.339/2],
 [0,0.339/2,0]]
 )

w = 0.3

B = np.array([
[0],
[-0.001017+2.217*w],
[-343.9/5+6.745*(w**2)],
[-0.2959/5+9.349*(w**3)],
[-1168+11.58*(w**4)],
[-1.663/2+8.68*(w**5)]])

[KD, KI, KP] = np.linalg.lstsq(A,B, rcond=None)[0]

# print(KD, KP, KI)
# print(w)

# Outro role 

w = 3

# Pra achar KI, na eq de s^2
KI = np.linalg.lstsq(np.array([[339/2000]]), np.array([[8.68*(w**5) - 1663/2000]]), rcond=None)[0]
print(KI)