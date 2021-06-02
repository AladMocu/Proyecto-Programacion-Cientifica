'''
=============================================================
                ECUACIONES DEL MODELO SEIR
=============================================================

    ds(t)/dt = -a_e s(t)e(t) - a_i s(t) i(t) + 𝛾 r(t)               (F1)
    de(t)/dt = a_e s(t)e(t) + a_i s(t)i(t) - ke(t) - pe(t)          (F2)
    di(t)/dt = ke(t) - 𝛽i(t) - 𝜇i(t)                                (F3)
    dr(t)/dt = 𝛽i(t) + 𝜌e(t) + 𝛾r(t)                                (F4)
    dp(t)/dt = 𝜇i(t)                                                (F5)

Con los siguientes parámetros:

param1: k: Tasa en la que aparecen los síntomas en los casos expuestos. 
param2: a_i: Factor de contagio entre las poblaciones infectadas y susceptibles.
param3: a_e: El factor de contagio entre las poblaciones susceptibles y las expuestas.
param4: y: Tasa de reinfección.
param5: b: La tasa de recuperación de los casos infectados.
param6: rho: La tasa de recuperación de los casos expuestos.
param7: mu: Tasa de mortalidad de los casos infectados.

Además, definimos s(t), e(t), i(t), r(t), p(t) de la siguiente manera:

s: s(t): Fracción de la población que es susceptible.
e: e(t): Fracción de la población que es expuesta.
i: i(t): Fracción de la población infectada.
r: r(t): Fracción de la población recuperada.
p: p(t): El número de individuos que muere debido a la enfermedad.
'''

import numpy as np
from scipy.integrate import odeint
import scipy.optimize as opt

iniciales = [0.8, 0.03, 0.03, 0.04, 0.1]
h = 0.1


# Definimos la función F1 ds(t)/dt:
def F1(s, e, i, r, a_e, a_i, y):
    return -a_e * s * e - a_i * s * i + y * r


# Definimos la función F2 de(t)/dt:
def F2(s, e, i, a_e, a_i, k, rho):
    return a_e * s * e + a_i * s * i - k * e - rho * e


# Definimos la función F3 di(t)/dt:
def F3(e, i, k, b, mu):
    return k * e - b * i - mu * i


# Definimos la función F4 dr(t)/dt:
def F4(i, e, r, b, rho, y):
    return b * i + rho * e + y * r


# Definimos la función F5 dp(t)/dt:
def F5(mu, i):
    return mu * i


'''
===========================================================================
                 SOLUCIÓN DEL SISTEMA DE ECUACIONES
===========================================================================

'''


# --------------------   FUNCIONES AUXILIARES    --------------------------

# EULER BACKWARD:
# Función cun sistemas de EDOs, que se usará posteriormente
# para resolver el sistema con fsolve.
# param1: aux: Matriz que almacena los vectores de solución de la iteración actual y[i] para
# s, e, i, r, p.
# param1: st1: Arreglo de soluciones de la fracción de la población susceptible, en la iteración anterior.
# param2: et1: Arreglo de soluciones de la fracción de la población expuesta, en la iteración anterior.
# param3: it1: Arreglo de soluciones de la fracción de la población infectada, en la iteración anterior.
# param4: rt1: Arreglo de soluciones de la fracción de la población recuperada, en la iteración anterior.
# param5: pt1: Arreglo de soluciones del número de individuos que muere debido a la enfermedad, en la iteración anterior.
# params: a_e, a_i, k, y, b, rho, mu.
# param6: h:

def FEulerBackRoot(aux, st1, et1, it1, rt1, pt1, a_e, a_i, k, y, b, rho, mu, h):
    return [st1 + h * F1(aux[0], aux[1], aux[2], aux[3], a_e, a_i, y) - aux[0],
            et1 + h * F2(aux[0], aux[1], aux[2], a_e, a_i, k, rho) - aux[1],
            it1 + h * F3(aux[1], aux[2], k, b, mu) - aux[2],
            rt1 + h * F4(aux[2], aux[3], b, rho, y) - aux[3],
            pt1 + h * F5(aux[2], mu) - aux[4]]

# EULER MODIFICADO:
# Función cun sistemas de EDOs, que se usará posteriormente
# para resolver el sistema con fsolve.
# param1: aux: Matriz que almacena los vectores de solución de la iteración actual y[i] para
# s, e, i, r, p.
# param1: st1: Arreglo de soluciones de la fracción de la población susceptible, en la iteración anterior.
# param2: et1: Arreglo de soluciones de la fracción de la población expuesta, en la iteración anterior.
# param3: it1: Arreglo de soluciones de la fracción de la población infectada, en la iteración anterior.
# param4: rt1: Arreglo de soluciones de la fracción de la población recuperada, en la iteración anterior.
# param5: pt1: Arreglo de soluciones del número de individuos que muere debido a la enfermedad, en la iteración anterior.
# params: a_e, a_i, k, y, b, rho, mu.
# param6: h:

def FEulerModRoot(aux, st1, et1, it1, rt1, pt1, a_e, a_i, k, y, b, rho, mu, h):
    return [st1 + (h / 2.0) * \
            (F1(st1, et1, it1, rt1, pt1, a_e, a_i, y) + \
             F1(aux[0], aux[1], aux[2], aux[3], aux[4], a_e, a_i, y)) - aux[0],
            et1 + (h / 2.0) * \
            (F2(st1, et1, it1, a_e, a_i, k, rho) + \
             F2(aux[0],aux[1], aux[2], a_e, a_i, k, rho)) - aux[1],
            it1 + (h / 2.0) * \
            (F3(et1, it1, k, b, mu) + \
             F3(aux[1], aux[2], k, b, mu)) - aux[2],
            rt1 + (h / 2.0) * \
            (F4(et1, it1, rt1, b, rho, y) + \
             F4(aux[1], aux[2], aux[3], b, rho, y)) - aux[3],
            pt1 + (h / 2.0) * \
            (F5(it1, u) + \
             F5(aux[2], yt2[3], mu)) - aux[4],
            ]

#-------------------------    MÉTODOS NUMÉRICOS   ---------------------------

#Método de condiciones iniciales:
def init_arr(time):
    s = np.zeros(len(time))
    e = np.zeros(len(time))
    i = np.zeros(len(time))
    r = np.zeros(len(time))
    p = np.zeros(len(time))
    #Valores iniciales teniendo en cuenta que se debe cumplir:
    #      s(t) + e(t) + i(t) + r(t) + p(t) = 1
    s[0] = iniciales[0]
    e[0] = iniciales[1]
    i[0] = iniciales[2]
    r[0] = iniciales[3]
    p[0] = iniciales[4]

    return s, e, i, r, p

#EULER FORWARD:

def euler_forward(params, time):
    a_e, a_i, k, y, b, rho, mu = params
    S_EulerFor, E_EulerFor, I_EulerFor, R_EulerFor, P_EulerFor = init_arr(time)
    h = abs(time[1] - time[0])

    for i in range(1, len(time)):
        S_EulerFor[iter] = S_EulerFor[i-1] + h * F1(S_EulerFor[i-1], E_EulerFor[i-1], I_EulerFor[i-1],
                                                   R_EulerFor[i-1], P_EulerFor[i-1], a_e, a_i, y)
        E_EulerFor[iter] = E_EulerFor[i-1] + h * F2(S_EulerFor[i-1], E_EulerFor[i-1], I_EulerFor[i-1],
                                                  a_e, a_i, k, rho)
        I_EulerFor[iter] = I_EulerFor[i-1] + h * F3(E_EulerFor[i-1], I_EulerFor[i-1], k, b, mu)
        R_EulerFor[iter] = R_EulerFor[i-1] + h * F4(I_EulerFor[i-1], E_EulerFor[i-1], R_EulerFor[i-1],
                                                    b, rho, y)
        P_EulerFor[iter] = P_EulerFor[i-1] + h * F5(I_EulerFor[i-1], mu)

    return S_EulerFor, E_EulerFor, I_EulerFor, R_EulerFor, P_EulerFor

#EULER BACKWARD
def euler_backward(params, time):
    a_e, a_i, k, y, b, rho, mu = params
    S_EulerBack, E_EulerBack, I_EulerBack, R_EulerBack, P_EulerBack = init_arr(time)

    for i in range(1, len(time)):
        SolBack = opt.fsolve(FEulerBackRoot, np.array([S_EulerBack[i-1], E_EulerBack[i-1], I_EulerBack[i-1],
                                                       R_EulerBack[i-1], P_EulerBack[i-1]]),
                             (S_EulerBack[i-1], E_EulerBack[i-1], I_EulerBack[i-1], R_EulerBack[i-1],
                             P_EulerBack[i-1], a_e, a_i, k, y, b, rho, mu))
        S_EulerBack[i] = SolBack[0]
        E_EulerBack[i] = SolBack[1]
        I_EulerBack[i] = SolBack[2]
        R_EulerBack[i] = SolBack[3]
        P_EulerBack[i] = SolBack[4]
    return S_EulerBack, E_EulerBack, I_EulerBack, R_EulerBack, P_EulerBack

#EULER MODIFICADO
def euler_modified(params, time):
    a_e, a_i, k, y, b, rho, mu = params
    S_EulerMod, E_EulerMod, I_EulerMod, R_EulerMod, P_EulerMod = init_arr(time)

    for i in range(1, len(time)):
        SolMod = opt.fsolve(FEulerModRoot, np.array([S_EulerMod[i - 1], E_EulerMod[i - 1], I_EulerMod[i - 1],
                                                       R_EulerMod[i - 1], P_EulerMod[i - 1]]),
                             (S_EulerMod[i - 1], E_EulerMod[i - 1], I_EulerMod[i - 1], R_EulerMod[i - 1],
                              P_EulerMod[i - 1], a_e, a_i, k, y, b, rho, mu))
        S_EulerMod[i] = SolMod[0]
        E_EulerMod[i] = SolMod[1]
        I_EulerMod[i] = SolMod[2]
        R_EulerMod[i] = SolMod[3]
        P_EulerMod[i] = SolMod[4]
    return S_EulerMod, E_EulerMod, I_EulerMod, R_EulerMod, P_EulerMod

# RK2
def runge_2(params, time):
    k, ai, ae, g, b, rho, u = params
    s, e, i, r, p = init_arr(time)

    h = abs(time[1] - time[0])
    for it in range(1, len(time)):
        ks1 = F1(s[it - 1], e[it - 1], i[it - 1], r[it - 1], ae, ai, g)
        ks2 = F1(s[it - 1] + ks1 * h, e[it - 1] + ks1 * h, i[it - 1] + ks1 * h, r[it - 1] + ks1 * h, ae, ai, g)
        s[it] = s[it - 1] + (h / 2) * (ks1 + ks2)

        ke1 = F2(s[it - 1], e[it - 1], i[it - 1], ae, ai, k, rho)
        ke2 = F2(s[it - 1] + ke1 * h, e[it - 1] + ke1 * h, i[it - 1] + ke1 * h, ae, ai, k, rho)
        e[it] = e[it - 1] + (h / 2) * (ke1 + ke2)

        ki1 = F3(e[it - 1], i[it - 1], k, b, u)
        ki2 = F3(e[it - 1] + ki1 * h, i[it - 1] + ki1 * h, k, b, u)
        i[it] = i[it - 1] + (h / 2) * (ki1 + ki2)

        kr1 = F4(e[it - 1], i[it - 1], r[it - 1], b, rho, g)
        kr2 = F4(e[it - 1] + kr1 * h, i[it - 1] + kr1 * h, r[it - 1] + kr1 * h, b, rho, g)
        r[it] = r[it - 1] + (h / 2) * (kr1 + kr2)

        kp1 = F5(i[it - 1], u)
        kp2 = F5(i[it - 1] + kp1 * h,u)

        p[it] = p[it - 1] + (h / 2) * (kp1 + kp2)

    return s, e, i, r, p


# RK4
def runge_4(params, time):
    k, ai, ae, g, b, rho, u = params
    s, e, i, r, p = init_arr(time)

    h = abs(time[1] - time[0])
    for it in range(1, len(time)):
        ks1 = F1(s[it - 1], e[it - 1], i[it - 1], r[it - 1], ae, ai, g)
        ks2 = F1(s[it - 1] + 0.5 * h * ks1, e[it - 1] + 0.5 * h * ks1, i[it - 1] + 0.5 * h * ks1,
                 r[it - 1] + 0.5 * h * ks1, ae, ai, g)
        ks3 = F1(s[it - 1] + 0.5 * h * ks2, e[it - 1] + 0.5 * h * ks2, i[it - 1] + 0.5 * h * ks2,
                 r[it - 1] + 0.5 * h * ks2, ae, ai, g)
        ks4 = F1(s[it - 1] + h * ks3, e[it - 1] + h * ks3, i[it - 1] + h * ks3, r[it - 1] + h * ks3, ae, ai, g)
        s[it] = s[it - 1] + (h / 6.0) * (ks1 + 2.0 * ks2 + 2.0 * ks3 + ks4)

        ke1 = F2(s[it - 1], e[it - 1], i[it - 1], ae, ai, k, rho)
        ke2 = F2(s[it - 1] + 0.5 * h * ke1, e[it - 1] + 0.5 * h * ke1, i[it - 1] + 0.5 * h * ke1, ae, ai, k, rho)
        ke3 = F2(s[it - 1] + 0.5 * h * ke2, e[it - 1] + 0.5 * h * ke2, i[it - 1] + 0.5 * h * ke2, ae, ai, k, rho)
        ke4 = F2(s[it - 1] + h * ke3, e[it - 1] + h * ke3, i[it - 1] + h * ke3, ae, ai, k, rho)
        e[it] = e[it - 1] + (h / 6.0) * (ke1 + 2.0 * ke2 + 2.0 * ke3 + ke4)

        ki1 = F3(e[it - 1], i[it - 1], k, b, u)
        ki2 = F3(e[it - 1] + 0.5 * h * ki1, i[it - 1] + 0.5 * h * ki1, k, b, u)
        ki3 = F3(e[it - 1] + 0.5 * h * ki2, i[it - 1] + 0.5 * h * ki2, k, b, u)
        ki4 = F3(e[it - 1] + h * ki3, i[it - 1] + h * ki3, k, b, u)
        i[it] = i[it - 1] + (h / 6.0) * (ki1 + 2.0 * ki2 + 2.0 * ki3 + ki4)

        kr1 = F4(e[it - 1], i[it - 1], r[it - 1], b, rho, g)
        kr2 = F4(e[it - 1] + 0.5 * h * kr1, i[it - 1] + 0.5 * h * kr1, r[it - 1] + 0.5 * h * kr1, b, rho, g)
        kr3 = F4(e[it - 1] + 0.5 * h * kr2, i[it - 1] + 0.5 * h * kr2, r[it - 1] + 0.5 * h * kr2, b, rho, g)
        kr4 = F4(e[it - 1] + h * kr3, i[it - 1] + h * kr3, r[it - 1] + h * kr3, b, rho, g)
        r[it] = r[it - 1] + (h / 6.0) * (kr1 + 2.0 * kr2 + 2.0 * kr3 + kr4)

        kp1 = F5(i[it - 1], u)
        kp2 = F5(i[it - 1] + 0.5 * h * kp1, u)
        kp3 = F5(i[it - 1] + 0.5 * h * kp2, u)
        kp4 = F5(i[it - 1] + h * kp3, u)
        p[it] = p[it - 1] + (h / 6.0) * (kp1 + 2.0 * kp2 + 2.0 * kp3 + kp4)

    return s, e, i, r, p


def aux_odeint(z, t, k, ai, ae, y, b, pp, u):
    s, e, i, r, p = z
    dsdt = -ae * s * e - ai * s * i + y * r
    dedt = ae * s * e + ai * s * i - k * e - pp * e
    didt = k * e - b * i - u * i
    drdt = b * i + p * e - y * r
    dpdt = u * i
    return [dsdt, dedt, didt, drdt, dpdt]


def odeint_s(params, range):

    z = odeint(aux_odeint, iniciales, range, args=tuple(params))
    return z[:, 0], z[:, 1], z[:, 2], z[:, 3], z[:, 4]


def solve(method, params, range):
    return {
        'Euler Forward': euler_forward(params, range),
        'Euler Backward': euler_backward(params, range),
        'Euler Modified': euler_modified(params, range),
        'Runge-Kutta 2': runge_2(params, range),
        'Runge-Kutta 4': runge_4(params, range),
        'odeint/ivp-solve': odeint_s(params, range),
    }[method]


if __name__ == "__main__":
    print("por favor ejecuta main.py")
