import numpy as np
import ctypes

#load the library
fun = ctypes.CDLL('./functions.so')
fun.F2sym.argtypes = [ctypes.c_float,ctypes.c_float,ctypes.c_float]
#define the result type of function
fun.F2sym.restype = ctypes.c_float

def pyF2sym(var):
    # To pass variables, wrap them as c variables 
    k = var[0]
    q = var[1]
    mq = var[2]
    return fun.F2sym(k,q,mq)/1e5 

def F2(var):
    k = var[0]
    q = var[1]
    mq = var[2]
    temp = (k**2*(7*k*mq + (3-10*mq**2)*q))/(14*q*(k**2 -2*k*mq*q + q**2))
    return temp


if __name__ =="__main__":
    u = [0.39188816, 0.87543167, 0.87426436]
#     u = [5,2,3]
    print( pyF2sym(u) )
    print( F2(u) )
