import numpy as np


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
