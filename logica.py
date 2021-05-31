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
                        CONDICIONES DEL SISTEMA
===========================================================================
'''





'''
===========================================================================
                 SOLUCIÓN DEL SISTEMA DE ECUACIONES
===========================================================================
'''
#Función para encontrar simultáneamente las raíces de las
#ecuaciones resultantes del método de Euler hacia atrás,
#para y1, y2, y3, y4, y5.
#param yt1:
#param yt2:
#param yt3:
#param yt4:
#param yt5:
# y_i = y_(i-1) + h * F(y_i)    →   0 = y_(i-1) + h * F(y_i) - y_i

def FEulerBackRoot(yt2, y1t1, y2t1, y3t1, y4t1, y5t1, h, k, a_i, a_e, g, b, p, u):
    return [y1t1 + h * F1(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, g) - yt2[0],
            y2t1 + h * F2(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], a_e, a_i, k, p) - yt2[1],
            y3t1 + h * F3(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], k, b, u) - yt2[2],
            y4t1 + h * F4(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4], b, p, g) - yt2[3],
            y5t1 + h * F5(yt2[0], yt2[1], yt2[2], yt2[3], yt2[4],  u) - yt2[4]]

#Función para encontrar simultáneamente las raíces de las
#ecuaciones resultantes del método de Euler modificado,
#para y1, y2, y3, y4, y5.
#param yt1:
#param yt2:
#param yt3:
#param yt4:
#param yt5:

def FEulerModRoot(yt2, y1t1, y2t1, h, k, a_i, a_e, g, b, p, u):
    retrurn [y1t1 + (h / 2.0) * ()]



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
