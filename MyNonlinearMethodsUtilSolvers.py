from math import *
import re
from scipy import integrate,interpolate,optimize
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
import numdifftools
import numpy as np
import matplotlib.pyplot as plt
import numbers
import decimal
import math
#TODO: include pyqtgraph, matplotlib for visualizations of methods


def bisectionExample(*,f=lambda x: 4*x**3 - 6*x**2 + 7*x - 2.3,a0=0,b0=1,tol=0.005,n=None):
    """
    Bisection (Bracketing) Method example
    
    f : #equation
    a0,b0: #NOTE: arbitrarily chosen by plot analysis
    tol: #epsilon tolerance, absolute_err_i <= tol
    n: #least set number of iterations. #NOTE: if precised (not None), tol is ignored.
    """
    assert (tol is None) ^ (n is None)
    
    if n is None:
        n = ceil(log2(abs(b0 - a0)/tol) - 1) #number of iterations, err <= |b0-a0|/2**n+1 <= tol

    mem = [[0,0] for _ in range(n+1)] # intervals [ai,bi]
    mem[0] = [a0,b0]

    xi = lambda i: sum(mem[i])/2 #xi+1 function for bisection method
    for i in range(n):
        x = xi(i)
        fx = f(x)
        fa = f(mem[i][0])
        print("{} : f({:11.7g}) = {:.16f}".format(i,x,fx))
        
        if (fa*fx) < 0.:
                a,b = mem[i][0],x
        else:
                a,b = x, mem[i][1]
        mem[i+1] = [a,b]

    xr = xi(n)
    print(f"approx  x={xr},\n\twhere {f(xr)} ~ 0",end="\n\r")
    if tol is not None:
        print(f"\t(epsilon of {tol}: {'' if abs(0 - f(xr)) <= tol else 'NOT'} respected)")

def falsePositionExample(*,f = lambda x: x**3.5 - 80,a0=2,b0=5,tol=0.8):
    """
    False Position (Bracketing) Method Example
    
    
    f : #equation
    a0,b0: #NOTE: arbitrarily chosen by plot analysis
    tol: #epsilon tolerance, absolute_err_i <= tol
    """
    
    i = 0 #number of CURRENT iterations
    mem = {} # intervals [ai,bi]
    mem[0] = [a0,b0]
    _E = lambda i: abs(mem[i][1]-mem[i][0])/2#abs error inequation of bracketing methods

    xc = lambda a,b,fa,fb: (a*fb - b*fa)/(fb - fa)  #xi+1 function for false position method
    xi = lambda i: xc(*mem[i],f(mem[i][0]),f(mem[i][1]))
    while _E(i) > tol:
           
        fa = f(mem[i][0])
        x = xi(i)
        fx = f(x)
        print("{} : f({:11.7g}) = {:.16f}".format(i,x,fx))
        if (fa*fx) < 0.:
                a,b = mem[i][0],x
        else:
                a,b = x, mem[i][1]
        mem[i+1] = [a,b]
        i += 1

    xr = xi(i)
    
    print(f"approx  x={xr},\n\twhere {f(xr)} ~ 0\n\t\
    (epsilon of {tol}: {'' if abs(0 - f(xr)) <= tol else 'NOT'} respected)"
    )

def fixedPointExample():
    """ 
    Fixed Point (Open) Method Example
    """
    #TODO demonstrate with 2+ different g(x) formula 
    #TODO and choose the fastest converging one
def newtonsExample(*,f = lambda x: 4*x**3 - 6*x**2 + 7*x - 2.3, tol = 0.005, x0 = 0, a0 = None , b0 = None):
    """
    Newton's (Open) Method Example
    
    f : #equation
    a0,b0: Interval in which root should be found. if precised, x0 is ignored and a matplotlib plot is launched in order to choose a x0.
    x0 : Initial guess. Ignored if a0 AND b0 precised.
    tol: #epsilon tolerance, absolute_err_i <= tol
    """
    
    if isinstance(a0,numbers.Number) and isinstance(b0,numbers.Number):
        DOM_SIZE = abs(b0 - a0)*(15/11)
        a0,b0 = min(b0,a0),max(b0,a0)
        VIEWDOM = np.linspace(a0 - DOM_SIZE*(2/11), 
                              b0 + DOM_SIZE*(2/11),
                              100)
        fv = np.vectorize(f)
        VIEWRANGE = fv(VIEWDOM)
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('zero')
        
        
        plt.plot(VIEWDOM,VIEWRANGE,'b')
        plt.axvspan(a0,b0,color='y',alpha=0.4)
        plt.show()
        
        x0 = None
        
        
        while x0 is None or not (a0 <= x0 <= b0):
            x0 = input("[Newton-Raphson method] Enter your starting guess x0: ")
            x0 = float(x0) if isinstance(x0,numbers.Number) else (None if re.match(r'^-?\d+(?:\.\d+)$',x0) is None else float(x0) )
    
    
    
    
    _x = sym.Symbol('_x')
    f =  lambda x: x*cos(x) - 2*x**2 + 3*x - 1#!!!!!!!!!!!
    
    df_dx = sym.lambdify(_x,
     sym.diff(f(_x),_x) 
    ) #first derivative of continuously differentiable function f
    
    
    mem = {}
    mem[0] = x0

    
    xi_plus_1 = lambda xi: xi - f(xi)/df_dx(xi)
    i = -1
    _E = lambda i: abs(mem[i]-mem[i+1])#abs error inequation of open methods
    
    while not ~i or _E(i) > tol:
        i += 1
        x = mem[i]
        fx = f(x)
        print("{} : f({:11.7g}) = {:.16f}".format(i,x,fx))
        
        mem[i+1] = xi_plus_1(x)
    
    xr = mem[i]
    print(f"approx  x={xr},\n\twhere {f(xr)} ~ 0\n\t\
    (epsilon of {tol}: {'' if abs(0 - f(xr)) <= tol else 'NOT'} respected)"
    )


def NewtonInterpPolyCoeff(x, y, num_of_coeffs_wanted, current_i=0):
    n = len(x)
    coefficients = [0] * (n - current_i)

    for i in range(current_i, n):
        divided_diff = y[i]
        for j in range(current_i, i):
            divided_diff = (divided_diff - coefficients[j - current_i]) / (x[i] - x[j])
        coefficients[i - current_i] = divided_diff

    return coefficients[:num_of_coeffs_wanted]

print(NewtonInterpPolyCoeff([1.0,2.0,3.0,4.0],[4.3,-2.6,-13.5,-70.4], 4))
class ComplexTrigInterp:

    def __init__(self,extraFloatPrec=False):
        self.__arginterp = 2
        self.trunc_precerror_maxbit = 128 if extraFloatPrec else 64 #TODO : Allocate the possibility of 32 bits for older architectures !
    pass
if __name__ == '__main__':
    #newtonsExample()
    eqsKwargs1 = {
        'f' : lambda x: sin(e**(3*x) + 4) - 0.5,
        'a0' : 0,
        'b0': .5,
        'tol' : .9e-10,
    }
    
    eqsKwargs2 = {
        'f' : lambda x: 4*cos(4*x)*sin(5*x) - 3,
        'a0' : 1.5,
        'b0': 2,
        'tol' : .9e-14,
    }
    
    eqsKwargs3 = {
        'f' : lambda x: x**3 * sin(3*x+4)-4,
        'a0' : 7,
        'b0': 8,
        'tol' : .9e-12,
    }
    #1.8293836019
    eqsKwargs4 = {
        'f' : lambda x: (x-2)**2 - log(x),
        'a0' : 1,
        'b0': 2,
        'tol' : .9e-9,
    }

    print("BISECTION")
    bisectionExample(**eqsKwargs1,n=None)
    input("Click to Proceed to next method")
    print("FALSE POINT")
    falsePositionExample(**eqsKwargs1)

    # eqsKwargs = {
    #     'f' : lambda x: x*cos(x) - 2*x**2 + 3*x - 1,
    #     'a0' : 1,
    #     'b0': 2,
    #     'tol' : .9e-9,
    # }
    # print("NEWTON")
    # newtonsExample(**eqsKwargs)