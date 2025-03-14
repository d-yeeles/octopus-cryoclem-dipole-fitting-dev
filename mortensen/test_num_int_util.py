"""
Python 3 Compatible Numerical Integration Utilities for PSF Models

Contains utility functions for numerical integration and other math operations
needed by both PSF models.

Updated to Python 3: February 25, 2025
"""

from scipy.special import jn, gamma, erfc, i1
from scipy.optimize import fmin_powell
from numpy import *
import numpy as np

###############################################################################
# Numerical Integration Utilities
###############################################################################

def trapzd(func, a, b, n):
    """
    Computes the n'th stage of refinement of an extended trapezoidal rule.

    With n=1 the function returns the crudest estimate of the integral of f
    in the interval [a,b]. Subsequent calls with n=2,3,... will improve the
    accuracy by adding 2^(n-2) additional interior point.
    """
    spacing = (b-a)/(2**n)
    xvals = arange(a+0.5*spacing, b, spacing)
    funcvals = func(xvals)
    funcsum = sum(funcvals)
    val = (b-a)*funcsum/(2**n)
    return val
            
def qtrap(func, a, b, eps=1.0e-10):
    """
    Integration is performed via the trapezoidal rule until a fractional
    accuracy of eps is obtained. Calls the routine trapzd.
    """
    jmax = 200
    val = oldval = 0.0
    
    for j in range(jmax):
        val = trapzd(func, a, b, j+1)
        if (j > 5):  # To avoid spurious early convergence
            if (fabs(val-oldval) < eps*fabs(oldval) or (val == 0.0 and oldval == 0.0)):
                return val
        oldval = val

    print("Too many steps in routine qtrap")
    return val

def polint(xa, ya, x):
    """
    Given arrays xa[0,...,n-1] and ya[0,...,n-1] and a value x, this routine
    returns a value y and an error estimate dy. If P(x) is the polynomial of
    degree N-1 such that P(xa[i])=ya[i] for i=0,...,n-1, then the returned
    value is y=P(x).
    """
    n = len(xa)
    c = zeros(n, float)
    d = zeros(n, float)

    ns = 0
    dif = fabs(x-xa[0])

    for i in range(n):
        dift = fabs(x-xa[i])
        if (dift < dif):  # Find the index of the closest table entry
            ns = i    
            dif = dift

        c[i] = ya[i]  # Initialize tableau of c's and d's
        d[i] = ya[i]

    y = ya[ns]  # This is the initial approximation to y.
    ns = ns-1

    for m in range(1, n):  # For each column in the tableau
        for i in range(n-m):
            ho = xa[i]-x
            hp = xa[i+m]-x

            w = c[i+1]-d[i]

            den = ho-hp
            if (den == 0.0):
                print("Error in routine polint")
                return (0, 0)

            den = w/den
            d[i] = hp*den
            c[i] = ho*den

        if (2*(ns+1) <= (n-m)):
            dy = c[ns+1]
        else:
            dy = d[ns]
            ns = ns-1

        y += dy
        
    return (y, dy)

def qromb(func, a, b, eps=1e-10):
    """
    Calculates the integral of the function f from a to b using 
    Romberg integration method.
    """
    jmax = 200
    jmaxp = jmax+1

    K = 5  # Number of points used in the extrapolation.

    s = zeros(jmax, float)
    s_t = zeros(K, float)
    h = zeros(jmaxp, float)
    h_t = zeros(K, float)

    h[0] = 1.0

    for j in range(jmax):
        s[j] = trapzd(func, a, b, j+1)
        if ((j+1) >= K):
            for i in range(K):
                h_t[i] = h[j+1-K+i]
                s_t[i] = s[j+1-K+i]

            (ss, dss) = polint(h_t, s_t, 0.0)
            if (fabs(dss) <= eps*fabs(ss)):
                return ss

        h[j+1] = 0.25*h[j]

    print("Too many steps in routine qromb!")
    return ss

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def fac(n):
    return gamma(n+1)

def erf(x):
    return 1-erfc(x)

def gauss(x, s):
    return exp(-x**2/(2*s**2))/sqrt(2*pi*s**2)
