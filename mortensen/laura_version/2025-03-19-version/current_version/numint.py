# Script containing functions to perform Romberg Integration on the real axis.
# The functions used as building blocks can be found in Numerical Recipes.

from numpy import arange, zeros, fabs

#from scipy.special import jn
#from simpson import Simpson

def trapzd(func,a,b,n):
    """
    Computes the n'th stage of refinement of an extended trapezoidal rule.

    With n=1 the function returns the crudest estimate of the integral of f
    in the interval [a,b]. Subsequent calls with n=2,3,... will improve the
    accuracy by adding 2^(n-2) additional interior point.

    See Numerical Recipes for a reference.

    Input:
    -------------------------------------------------------------------------
    func (real) : A function on the real axis returning its functional value.
    a,b (real)  : The limits of the integration.
    n (int)     : The order of refinement for the method.

    Output:
    -------------------------------------------------------------------------
    val (real)  : The value of the integral.
    """

    spacing=(b-a)/(2**n)
    xvals=arange(a+0.5*spacing,b,spacing)
    funcvals=func(xvals)
    funcsum=sum(funcvals)
    val=(b-a)*funcsum/(2**n)
    return val


def qtrap(func,a,b,eps=1.0e-10):
    """
    Integation is performed via the trapezoidal rule util a frational
    accuracy of eps is obtained. Calls the routine trapzd.

    See Numerical recipes for a reference.

    Input:
    -------------------------------------------------------------------------
    func (real) : A function on the real axis returning its functional value.
    a,b (real)  : The limits of the integration.
    eps (real)  : The fractional accuracy to be obtained.

    Output:
    -------------------------------------------------------------------------
    val (real)  : The value of the integral.

    """

    jmax=200

    val=oldval=0.0

    for j in range(jmax):
        val=trapzd(func,a,b,j+1)
        if (j>5): # To avoid spurious early convergence
            if (fabs(val-oldval)<eps*fabs(oldval) or (val==0.0 and oldval==0.0)):
                return val
        oldval=val

    print("Too many steps in routine qtrap")
    return



def polint(xa,ya,x):
    """
    Given arrays xa[0,...,n-1] and ya[0,...,n-1] and a value x, this routine
    returns a value y and an error estimate dy. If P(x) is the polynomial of
    degree N-1 such that P(xa[i])=ya[i] for i=0,...,n-1, then the returned
    value is y=P(x).

    Input:
    -------------------------------------------------------------
    xa (array)  : Vector of length n containing x-values.
    ya (array)  : Vector of length n containing y-values.
    x (float)   : x values of which we need an extrapolated value.

    Output:
    -------------------------------------------------------------
    y (float)   : The extrapolated value at x.
    dy (float)  : An error estimate.
    """

    n=len(xa)
    c=zeros(n,float)
    d=zeros(n,float)

    ns=0
    dif=fabs(x-xa[0])

    for i in range(n):
        dift=fabs(x-xa[i])
        if (dift<dif):  # Find the index of the closest table entry
            ns=i
            dif=dift

        c[i]=ya[i] # Initialize tableau of c's and d's
        d[i]=ya[i]

    y=ya[ns] # This is the initial approximation to y.
    ns=ns-1

    for m in range(1,n): # For each collumn in the tableau
        for i in range(n-m):
            ho=xa[i]-x
            hp=xa[i+m]-x

            w=c[i+1]-d[i]

            den=ho-hp
            if (den==0.0):
                print("Error in routine polint")
                return

            den=w/den
            d[i]=hp*den
            c[i]=ho*den

        if (2*(ns+1)<=(n-m)):
            dy=c[ns+1]
        else:
            dy=d[ns]
            ns=ns-1

        y+=dy

    return (y,dy)

##xa=arange(2.0,6.5,0.7)
##ya=sin(xa)
##x=7.0
##(p,dp)=polint(xa,ya,x)
##a=arange(0,10,0.001)
##b=sin(a)
##plot(xa,ya,'bo',[x],[p],'go',a,b,'b-')


def qromb(func,a,b,eps=1e-10):
    """
    Calculates the integral of the function f from a to b. Integration is
    performed via the Romberg integration method of order 2K where K=2 is
    the Simpson's rule. The routine trapzd is called.

    eps is the fractional accuracy desired as determined by the extrapolation
    error estimate. jmax limits the total number of desired steps. K is the
    number of points used in the extrapolation. The routine polint is called.

    Input:
    -------------------------------------------------------------------------
    func (float) : A function on the real axis returning its functional value.
    a,b (float)  : The limits of the integration.
    eps (float)  : The fractional accuracy to be obtained.

    Output:
    -------------------------------------------------------------------------
    val (float)  : The value of the integral.

    """

    jmax=200
    jmaxp=jmax+1

    K=5 # Number of points used in the extrapolation.

    s=zeros(jmax,float)
    s_t=zeros(K,float)
    h=zeros(jmaxp,float)
    h_t=zeros(K,float)

    h[0]=1.0

    for j in range(jmax):
        s[j]=trapzd(func,a,b,j+1)
        if ((j+1)>=K):
            for i in range(K):
                h_t[i]=h[j+1-K+i]
                s_t[i]=s[j+1-K+i]

            (ss,dss)=polint(h_t,s_t,0.0)
            if (fabs(dss)<=eps*fabs(ss)):
                return ss

        h[j+1]=0.25*h[j]

    print("Too many steps in routine qromb!")
    return


##def test(x):
##    value=1.0/x
##    return value
##
##def test2(x):
##    value=x**4*log(x+sqrt(x**2+1.0))
##    return value
##
##ep=2.0
##
##print "log2 =",log(ep)
##print "qromb=",qromb(test,1.0,ep,1.0e-10)
##print "qtrap=",qtrap(test,1.0,ep,1.0e-10)
##print "simp =", Simpson(test,1.0,ep,m=1000)
##
##def lin(x):
##    val=2.0*x
##    return val
##
##def bessel_sq(x):
##    val=(jn(1,x))**2
##    return val

