import sympy as sym
import numpy as np
import types

s,c,alpha, l = sym.symbols("s,c,alpha,l")

def exp_omega_mat(vec):
    mat = sym.Matrix([[0,-vec[2], vec[1]],
                      [vec[2], 0, -vec[0]],
                      [-vec[1], vec[0], 0]])
    I = sym.eye(3)
    
    exp_mat = I + mat*s + (mat**2)*(1-c)
    
    vec_mat = sym.Matrix([[vec[0], vec[1], vec[2]]]).transpose()
    ww_mat = vec_mat*vec_mat.transpose()

    
    k = (I - exp_mat)*mat + alpha*ww_mat
    return k

def exp_twsit(twist):
    w = twist[3:]
    v = twist[:3]
    I = sym.eye(3)
    
    rot = exp_omega_mat(w)
    cross_product = (sym.Matrix(w).cross(sym.Matrix(v)))
    
    trans = (I - rot)*cross_product
    print("*"*50)
    print("rot :", rot)
    print(" trans :", trans)
    print("*"*50)
    
    
    

k  = sym.Matrix([[1, 2, alpha], [alpha,5,6], [7,2,6]])
twist1 = [0,0,0,0,0,1]
twist2= [l,0,0,0,0,1]
twist3 = [l+l,0,0,0,0,1]
twist4 = [0,0,1,0,0,0]

exp_twsit(twist1)
exp_twsit(twist2)
exp_twsit(twist3)
exp_twsit(twist4)


