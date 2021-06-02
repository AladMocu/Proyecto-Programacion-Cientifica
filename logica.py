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
param4: g: Tasa de reinfección.
param5: b: La tasa de recuperación de los casos infectados.
param6: p: La tasa de recuperación de los casos expuestos.
param7: u: Tasa de mortalidad de los casos infectados.

Además, definimos s(t), e(t), i(t), r(t), p(t) de la siguiente manera:

y1: s(t): Fracción de la población que es susceptible.
y2: e(t): Fracción de la población que es expuesta.
y3: i(t): Fracción de la población infectada.
y4: r(t): Fracción de la población recuperada.
y5: p(t): El número de individuos que muere debido a la enfermedad.
'''

import numpy as np
from scipy.integrate import odeint

iniciales = [0.8, 0.03, 0.03, 0.04, 0.1]
h = 0.1


# Definimos la función F1 ds(t)/dt:
def F1(s, e, i, r, a_e, a_i, y):
    return -a_e * s * e - a_i * s * i + y * r


# Definimos la función F2 de(t)/dt:
def F2(s, e, i, a_e, a_i, k, rho):
    return a_e * s * e + a_i * s * i - k * e - rho * e


# Definimos la función F3 di(t)/dt:
def F3(e, i, k, b, u):
    return k * e - b * i - u * i


# Definimos la función F4 dr(t)/dt:
def F4(i, e, r, b, rho, y):
    return b * i + rho * e + y * r


# Definimos la función F5 dp(t)/dt:
def F5(mu, u):
    return u * mu


'''
===========================================================================
                 SOLUCIÓN DEL SISTEMA DE ECUACIONES
===========================================================================

'''


# -----------------------FUNCIONES AUXILIARES-----------------------------

# EULER BACKWARD:
def FEulerBackRoot(params, time):
    a_e, a_i, k, g, b, p, u = params
    yt2, y1t1, y2t1, y3t1, y4t1, y5t1 = init_arr(time)
    h = abs(time[1] - time[0])
    return [y1t1 + h * F1(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, g) - yt2[0],
            y2t1 + h * F1(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, p) - yt2[1],
            y3t1 + h * F1(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], k, b, u) - yt2[2],
            y4t1 + h * F1(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], b, p, g) - yt2[3],
            y5t1 + h * F1(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], u) - yt2[4]]


# EULER MODIFICADO:
def FEulerModRoot(params, time):
    a_e, a_i, k, g, b, p, u = params
    yt2, y1t1, y2t1, y3t1, y4t1, y5t1 = init_arr(time)
    h = abs(time[1] - time[0])

    return [y1t1 + (h / 2.0) *
            (F1(y1t1, y2t1, y3t1, y4t1, y5t1, a_e, a_i, k, g, b, p, u) +
             F1(yt2[0],yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, g, b, p, u)) - yt2[0],
            y2t1 + (h / 2.0) *
            (F2(y1t1, y2t1, y3t1, y4t1, y5t1, a_e, a_i, k, g, b, p, u) +
             F2(yt2[0],yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, g, b, p, u)) - yt2[1],
            y3t1 + (h / 2.0) *
            (F3(y1t1, y2t1, y3t1, y4t1, y5t1, a_e, a_i, k, g, b, p, u) +
             F3(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, g, b, p, u)) - yt2[2],
            y4t1 + (h / 2.0) *
            (F4(y1t1, y2t1, y3t1, y4t1, y5t1, a_e, a_i, k, g, b, p, u) +
             F4(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, g, b, p, u)) - yt2[3],
            y5t1 + (h / 2.0) *
            (F5(y1t1, y2t1, y3t1, y4t1, y5t1, a_e, a_i, k, g, b, p, u) +
             F5(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, g, b, p, u)) - yt2[4],
            ]

#------------------------------Métodos numéricos-------------------------------
#Método de condiciones iniciales:
def init_arr(time):
    s = np.zeros(len(time))
    e = np.zeros(len(time))
    i = np.zeros(len(time))
    r = np.zeros(len(time))
    p = np.zeros(len(time))
    # Falta definir los valores iniciales
    s[0] = iniciales[0]
    e[0] = iniciales[1]
    i[0] = iniciales[2]
    r[0] = iniciales[3]
    p[0] = iniciales[4]

    return s, e, i, r, p

#EULER FORWARD:

#EULER BACKWARDS:

#EULER MODIFICADO:

def euler_forward(params, range):
    k, ai, ae, y, b, p, u = params
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30 / (params[0] + np.exp(-0.5 * range))
    e = 20 / (params[1] + np.exp(-0.5 * range))
    i = 10 / (params[2] + np.exp(-0.5 * range))
    r = 90 / (params[3] + np.exp(-0.5 * range))
    p = 80 / (params[4] + np.exp(-0.5 * range))
    return s, e, i, r, p


def euler_backward(params, range):
    k, ai, ae, y, b, p, u = params
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30 / (params[0] + np.exp(-0.5 * range))
    e = 20 / (params[1] + np.exp(-0.5 * range))
    i = 10 / (params[2] + np.exp(-0.5 * range))
    r = 90 / (params[3] + np.exp(-0.5 * range))
    p = 80 / (params[4] + np.exp(-0.5 * range))
    return s, e, i, r, p


def euler_modified(params, range):
    k, ai, ae, y, b, p, u = params
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30 / (params[0] + np.exp(-0.5 * range))
    e = 20 / (params[1] + np.exp(-0.5 * range))
    i = 10 / (params[2] + np.exp(-0.5 * range))
    r = 90 / (params[3] + np.exp(-0.5 * range))
    p = 80 / (params[4] + np.exp(-0.5 * range))
    return s, e, i, r, p


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
        kp2 = F5(i[it - 1] + kp1 * h, u)

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
