import numpy as np

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

#Definimos la función F1:
def F1(y1, y2, y3, y4, a_e, a_i, g):
    return -a_e * y1 * y2 - a_i * y1 * y3 + g * y4

#Definimos la función F2:
def F2(y1, y2, y3, a_e, a_i, k, p):
    return a_e * y1 *y2 + a_i * y1 * y3 - k * y2 - p * y2

#Definimos la función F3:
def F3(y2, y3, k, b, u):
    return k * y2 - b * y3 - u * y3

#Definimos la función F4:
def F4(y3, y2, y4, b, p, g):
    return b * y3 + p * y2 + g * y4

#Definimos la función F5:
def F5(y3, u):
    return u * y3


'''
===========================================================================
                 SOLUCIÓN DEL SISTEMA DE ECUACIONES
===========================================================================


'''

def demo_method(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def euler_forward(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def euler_backward(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def euler_modified(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def runge_2(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def runge_4(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def odeint(params, range):
    s = []
    e = []
    i = []
    r = []
    p = []

    s = 30/(params[0]+np.exp(-0.5*range))
    e = 20/(params[1]+np.exp(-0.5*range))
    i = 10/(params[2]+np.exp(-0.5*range))
    r = 90/(params[3]+np.exp(-0.5*range))
    p = 80/(params[4]+np.exp(-0.5*range))
    return s, e, i, r, p


def solve(method, params, range):
    return {
        'Euler Forward': euler_forward(params, range),
        'Euler Backward': euler_backward(params, range),
        'Euler Modified': euler_modified(params, range),
        'Runge-Kutta 2': runge_2(params, range),
        'Runge-Kutta 4': runge_4(params, range),
        'odeint/ivp-solve': odeint(params, range),
    }[method]


if __name__ == "__main__":
    print("por favor ejecuta main.py")
