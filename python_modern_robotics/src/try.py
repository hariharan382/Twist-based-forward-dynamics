from scipy.optimize import minimize
import numpy as np
from scipy.integrate import odeint

m = 1220
k = 35600
g = 17.5
a = 450000
z0 = np.array([-.5, 0])

def d(z, t, m, k, g, a, b):
    return np.array([-1/m * (b*z[0] + k*z[1] + a*z[1]**3 + m*g), z[0]])

def func(b, z0, *args):
    _, x = odeint(d, z0, t, args=args+(b,)).T
    return -x.max()  # minimize negative max

cons = [{'type': 'ineq', 'fun': lambda b: b - 1000, 'jac': lambda b: 1},   # b > 1000
        {'type': 'ineq', 'fun': lambda b: 10000 - b, 'jac': lambda b: -1}, # b < 10000
        {'type': 'ineq', 'fun': lambda b: func(b, z0, m, k, g, a)}] # func(b) > 0 means x < 0

b0 = 10000
b_min = minimize(func, b0, args=(z0, m, k, g, a), constraints=cons)