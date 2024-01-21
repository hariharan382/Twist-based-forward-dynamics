
import numpy as np
from sympy import symbols
import sympy as sym
from sympy.abc import x
from sympy import Function, dsolve, Eq, pprint, diff,sin, cos,pi
from pprint import pprint
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#init variables  
J = 1
m=1
g =10
kp =1/4
kd = 1
q0 = pi/8
qdot0 = 0
qd = pi/2
l=1


def ode_func(q, t, J,m,l,kp,kd,g,qd):
    dqdt = [q[1], (1/J)*(kp*(qd- q) + kd(q[1]) -m*g*l*sin(q[0]) + m*g*l*sin(qd))]
    
    return  dqdt



# t =  symbols("t")
# q = sym.Function("q", real=True)(t)
# dqdt = q.diff(t)
# dq2dt2 = q.diff(t,2)
# diff_eqn = sym.Eq(J*dq2dt2+m*g*l*sin(q)-kp*(qd-q) + kd*(dqdt) - m*g*l*sin(qd),0)
# ode_func = lambda q,t :[q[1], (1/J)*(kp*(qd- q) + kd(q[1]) -m*g*l*sin(q[0]) + m*g*l*sin(qd))]
t_points = np.linspace(0,10,50)
solution = odeint(ode_func, [q0,qdot0], t_points, args = (J,m,g,l,kd,kp,qd))
# eqn = sym.lambdify((t), diff_eqn)
# print("eqn :", eqn)
# solution =  scipy.integrate.solve_ivp(eqn, (0,10), [pi/8, 0], t_eval = t_eval)
# ics= {q.subs(t,0):q0, q.diff(t).subs(t,0):qdot0}
# result = sym.dsolve(diff_eqn,ics = ics)
pprint(solution)





