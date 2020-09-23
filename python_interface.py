# from ctypes import cdll,c_float,byref
import ctypes

#load the library
fun = ctypes.CDLL('./functions.so')
fun.F2sym.argtypes = [ctypes.c_float,ctypes.c_float,ctypes.c_float]
#define the result type of function
fun.F2sym.restype = ctypes.c_float


#Run the function
print( fun.F2sym(5.,2.,3.) )
