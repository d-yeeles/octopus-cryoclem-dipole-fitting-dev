# Supplementary Software 
#
# for 
#
# Optimized measurements of separations and angles between intra-molecular fluorescent markers
#
# by
#
# Kim I. Mortensen, Jongmin Sung, Henrik Flyvbjerg, and James A. Spudich
#
# Python (python.org) code for the calculation of the theoretical PSF
# for fixed dipoles imaged close to the focus, maximum likelihood
# estimator for the position and orientation of the probe,
# and the calculation of the covariance matrix for the estimated
# parameters.
#
# Uses the modules matplotlib, scipy, numpy
#
# September 8, 2015
#
# Kim I. Mortensen

from scipy.special import jn,gamma,erfc,i1
from scipy.optimize import fmin_powell, minimize
from pylab import *
import numpy

# Script containing functions to perform Romberg Integration on the real axis.
# The functions used as building blocks can be found in Numerical Recipes.

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


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def fac(n):
    return gamma(n+1)

def erf(x):
    return 1-erfc(x)

def gauss(x,s):
    return exp(-x**2/(2*s**2))/sqrt(2*pi*s**2)


class dipdistr:

    """

    Calculates the theoretical point spread function (PSF) for fixed dipoles.
    
    The PSF is the distribution of photons in the image plane from
    a fixed dipole emitter close to the coverslip surface,
    when imaged by a large aperture objective close to the design focal plane.

    Input:
    -------------------------------------------------------------------------
    
    wl (float) : Wavelength [nm] of emission peak in buffer.
    NA (float) : Numerical aperture of the objective
    n (float) : Refractive index of glass/immersion oil
    n0 (float) : Refractive index of buffer (water)
    M (float) : Magnification of the objective
    
    Functions:
    -------------------------------------------------------------------------
    PSF_exact (float) : Takes the coordinates in the image, the angles of the
    probe, and the distance from the focus as arguments. 
    Returns the value of the exact PSF at that point.

    PSF_approx (float): Takes the coordinates in the image, the angles of the
    probe, and the distance from the focus as arguments. 
    Returns the value of the approximated PSF at that point.

    """

    def __init__(self,wl,n,n0,M,NA):

        self.wl=wl          #Wavelength in nm (emission peak in buffer)
        self.NA=NA          #Numerical aperture
        
        self.M=M            #Magnification
        self.n=n            #Refractive index of glass/oil
        self.np=1.0         #Refractive index of air in lab (vacuum)
        self.n0=n0          #Refractive index of sample medium (water)

        self.kp=2*pi/wl     #Wavevector amplitude in air (vacuum)
        self.k0=n0*self.kp  #Wavevector amplitude in buffer
        self.k=n*self.kp    #Wavevector amplitude in glass

        self.etapmed=n0/M   #Integration limits
        self.etapmax=NA/M

        # Calculate or load normalization constants
        try:
            self.norm=numpy.loadtxt('dipoletablenorm.dat')

        except:
             self.deltaz=-1.0
             self.norm=self.Normalization()
             numpy.savetxt('dipoletablenorm.dat',self.norm)

        # Tabulate the PSF approximations at this range of focal values 
        self.focusvals=arange(-150,150,1.0)
        self.nfoci=len(self.focusvals)
        nf=self.nfoci

        self.b0subreal_norm=zeros(nf)
        self.b0subreal_mean=zeros(nf)
        self.b0subreal_var=zeros(nf)
        self.b0subimag_norm=zeros(nf)
        self.b0subimag_mean=zeros(nf)
        self.b0subimag_var=zeros(nf)
        self.b0real_norm=zeros(nf)
        self.b0real_mean=zeros(nf)
        self.b0real_var=zeros(nf)
        self.b0imag_norm=zeros(nf)
        self.b0imag_mean=zeros(nf)
        self.b0imag_var=zeros(nf)

        self.b1subreal_norm=zeros(nf)
        self.b1subreal_mean=zeros(nf)
        self.b1subreal_var=zeros(nf)
        self.b1subimag_norm=zeros(nf)
        self.b1subimag_mean=zeros(nf)
        self.b1subimag_var=zeros(nf)
        self.b1real_norm=zeros(nf)
        self.b1real_mean=zeros(nf)
        self.b1real_var=zeros(nf)
        self.b1imag_norm=zeros(nf)
        self.b1imag_mean=zeros(nf)
        self.b1imag_var=zeros(nf)

        self.b2subreal_norm=zeros(nf)
        self.b2subreal_mean=zeros(nf)
        self.b2subreal_var=zeros(nf)
        self.b2subimag_norm=zeros(nf)
        self.b2subimag_mean=zeros(nf)
        self.b2subimag_var=zeros(nf)
        self.b2real_norm=zeros(nf)
        self.b2real_mean=zeros(nf)
        self.b2real_var=zeros(nf)
        self.b2imag_norm=zeros(nf)
        self.b2imag_mean=zeros(nf)
        self.b2imag_var=zeros(nf)
        
        
        for nf in range(self.nfoci):
            self.deltaz=self.focusvals[nf]       

        # Calculate the zeroth, first, and second moments of eta' used in the
        # cumulant approximation of the true PSF. Do this separately for the
        # sub- and supercritical regions. This calculation is only needed when
        # the class is initiated.
        
            def Integrand0_subreal(etap):
                    eta=self.Eta(etap)
                    integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))\
                               *cos(self.k*self.deltaz*cos(eta))
                    return integrand

            self.b0subreal_norm[nf]=qromb(Integrand0_subreal,0.0,self.etapmed,1e-4)
            self.b0subreal_mean[nf]=qromb(lambda etap: Integrand0_subreal(etap)*etap,0.0,self.etapmed,1e-4)\
                             /self.b0subreal_norm[nf]
            self.b0subreal_var[nf]=qromb(lambda etap: Integrand0_subreal(etap)*etap**2,0.0,self.etapmed,1e-4)\
                            /self.b0subreal_norm[nf]-self.b0subreal_mean[nf]**2
           
            def Integrand0_subimag(etap):
                    eta=self.Eta(etap)
                    integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))\
                               *sin(self.k*self.deltaz*cos(eta))
                    return integrand

            if abs(self.deltaz)>=1.0:
                self.b0subimag_norm[nf]=qromb(Integrand0_subimag,0.0,self.etapmed,1e-4)
                self.b0subimag_mean[nf]=qromb(lambda etap: Integrand0_subimag(etap)*etap,0.0,self.etapmed,1e-4)\
                                 /self.b0subimag_norm[nf]
                self.b0subimag_var[nf]=qromb(lambda etap: Integrand0_subimag(etap)*etap**2,0.0,self.etapmed,1e-4)\
                                /self.b0subimag_norm[nf]-self.b0subimag_mean[nf]**2

            def Integrand0_real(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           -etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b0real_norm[nf]=qromb(Integrand0_real,self.etapmed,self.etapmax,1e-4)
            self.b0real_mean[nf]=qromb(lambda etap: Integrand0_real(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                              /self.b0real_norm[nf]
            self.b0real_var[nf]=qromb(lambda etap: Integrand0_real(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                             /self.b0real_norm[nf]-self.b0real_mean[nf]**2

            def Integrand0_imag(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           +etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b0imag_norm[nf]=qromb(Integrand0_imag,self.etapmed,self.etapmax,1e-4)
            self.b0imag_mean[nf]=qromb(lambda etap: Integrand0_imag(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                              /self.b0imag_norm[nf]
            self.b0imag_var[nf]=qromb(lambda etap: Integrand0_imag(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                             /self.b0imag_norm[nf]-self.b0imag_mean[nf]**2

            def Integrand1_subreal(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*self.Epperp(etap)\
                           *cos(self.k*self.deltaz*cos(eta))
                return integrand

            self.b1subreal_norm[nf]=qromb(Integrand1_subreal,0.0,self.etapmed,1e-4)
            self.b1subreal_mean[nf]=qromb(lambda etap: Integrand1_subreal(etap)*etap,0.0,self.etapmed,1e-4)\
                             /self.b1subreal_norm[nf]
            self.b1subreal_var[nf]=qromb(lambda etap: Integrand1_subreal(etap)*etap**2,0.0,self.etapmed,1e-4)\
                            /self.b1subreal_norm[nf]-self.b1subreal_mean[nf]**2

            def Integrand1_subimag(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*self.Epperp(etap)\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            if abs(self.deltaz)>=1.0:
                self.b1subimag_norm[nf]=qromb(Integrand1_subimag,0.0,self.etapmed,1e-4)
                self.b1subimag_mean[nf]=qromb(lambda etap: Integrand1_subimag(etap)*etap,0.0,self.etapmed,1e-4)\
                                 /self.b1subimag_norm[nf]
                self.b1subimag_var[nf]=qromb(lambda etap: Integrand1_subimag(etap)*etap**2,0.0,self.etapmed,1e-4)\
                                /self.b1subimag_norm[nf]-self.b1subimag_mean[nf]**2

            def Integrand1_real(etap):
                eta=self.Eta(etap)
                integrand=integrand=etap/sqrt(cos(eta))*self.sc5(etap)\
                           *cos(self.k*self.deltaz*cos(eta))\
                           -etap/sqrt(cos(eta))*self.sc6(etap)\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b1real_norm[nf]=qromb(Integrand1_real,self.etapmed,self.etapmax,1e-4)
            self.b1real_mean[nf]=qromb(lambda etap: Integrand1_real(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                              /self.b1real_norm[nf]
            self.b1real_var[nf]=qromb(lambda etap: Integrand1_real(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                             /self.b1real_norm[nf]-self.b1real_mean[nf]**2

            def Integrand1_imag(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*self.sc6(etap)\
                           *cos(self.k*self.deltaz*cos(eta))\
                           +etap/sqrt(cos(eta))*self.sc5(etap)\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b1imag_norm[nf]=qromb(Integrand1_imag,self.etapmed,self.etapmax,1e-4)
            self.b1imag_mean[nf]=qromb(lambda etap: Integrand1_imag(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                              /self.b1imag_norm[nf]
            self.b1imag_var[nf]=qromb(lambda etap: Integrand1_imag(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                             /self.b1imag_norm[nf]-self.b1imag_mean[nf]**2


            def Integrand2_subreal(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))\
                           *cos(self.k*self.deltaz*cos(eta))
                return integrand

            self.b2subreal_norm[nf]=qromb(Integrand2_subreal,0.0,self.etapmed,1e-4)
            self.b2subreal_mean[nf]=qromb(lambda etap: Integrand2_subreal(etap)*etap,0.0,self.etapmed,1e-4)\
                             /self.b2subreal_norm[nf]
            self.b2subreal_var[nf]=qromb(lambda etap: Integrand2_subreal(etap)*etap**2,0.0,self.etapmed,1e-4)\
                            /self.b2subreal_norm[nf]-self.b2subreal_mean[nf]**2


            def Integrand2_subimag(etap):
                eta=self.Eta(etap)
                integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            if abs(self.deltaz)>=1.0:
                self.b2subimag_norm[nf]=qromb(Integrand2_subimag,0.0,self.etapmed,1e-4)
                self.b2subimag_mean[nf]=qromb(lambda etap: Integrand2_subimag(etap)*etap,0.0,self.etapmed,1e-4)\
                                 /self.b2subimag_norm[nf]
                self.b2subimag_var[nf]=qromb(lambda etap: Integrand2_subimag(etap)*etap**2,0.0,self.etapmed,1e-4)\
                                /self.b2subimag_norm[nf]-self.b2subimag_mean[nf]**2


            def Integrand2_real(etap):
                eta=self.Eta(etap)
                integrand=integrand=etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           -etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b2real_norm[nf]=qromb(Integrand2_real,self.etapmed,self.etapmax,1e-4)
            self.b2real_mean[nf]=qromb(lambda etap: Integrand2_real(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                              /self.b2real_norm[nf]
            self.b2real_var[nf]=qromb(lambda etap: Integrand2_real(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                             /self.b2real_norm[nf]-self.b2real_mean[nf]**2

            def Integrand2_imag(etap):
                eta=self.Eta(etap)
                integrand=integrand=etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           +etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b2imag_norm[nf]=qromb(Integrand2_imag,self.etapmed,self.etapmax,1e-4)
            self.b2imag_mean[nf]=qromb(lambda etap: Integrand2_imag(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                              /self.b2imag_norm[nf]
            self.b2imag_var[nf]=qromb(lambda etap: Integrand2_imag(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                             /self.b2imag_norm[nf]-self.b2imag_mean[nf]**2


    # Change and retrieve the focus

    def SetFocus(self,deltaz):
        self.deltaz=deltaz

    def GetFocus(self):
        return self.deltaz
 
    # Change and retrieve the radial position value
    
    def SetRho(self,rho):
        self.rho=rho

    def GetRho(self):
        return self.rho

    # Connect eta and eta0 to etap via Snell's law

    def Eta(self,etap):
        M=self.M
        np=self.np
        n=self.n
        eta=arcsin(M*np/n*etap)
        return eta

    def Eta0(self,etap):
        M=self.M
        np=self.np
        n0=self.n0
        n=self.n
        eta=self.Eta(etap)
        eta0=arccos(sqrt(fabs(1.0-(n/n0*sin(eta))**2)))
        return eta0
    
    # The Fresnel transmission coefficients
    
    def Ts(self,etap):
        eta=self.Eta(etap)
        eta0=self.Eta0(etap)
        ts=2*cos(eta0)*sin(eta)/sin(eta0+eta)
        return ts

    def Tp(self,etap):
        eta=self.Eta(etap)
        eta0=self.Eta0(etap)
        ts=self.Ts(etap)
        tp=ts/cos(eta0-eta)
        return tp

    # z-components of wavevectors

    def W(self,etap):
        eta=self.Eta(etap)
        k=self.k
        w=cos(eta)*k
        return w

    def W0(self,etap):
        eta0=self.Eta0(etap)
        k0=self.k0
        w0=cos(eta0)*k0
        return w0

    # Support functions

    def com(self,etap):
        n=self.n
        n0=self.n0
        eta=self.Eta(etap)
        value=sqrt(fabs(1-(n/n0*sin(eta))**2))
        return value

    def gamma(self,etap):
        n=self.n
        n0=self.n0
        eta=self.Eta(etap)
        c=self.com(etap)
        value=(n/n0)*cos(eta)/c
        return value

    def delta(self,etap):
        n=self.n
        n0=self.n0
        eta=self.Eta(etap)
        c=self.com(etap)
        value=(n0/n)*cos(eta)/c
        return value

    def epsilon(self,etap):
        n=self.n
        n0=self.n0
        eta=self.Eta(etap)
        c=self.com(etap)
        value=(n/n0)*c/cos(eta)
        return value

    # Integrands for super-critical angles

    def sc1(self,etap):
        n,n0=self.n,self.n0
        k,k0,kp=self.k,self.k0,self.kp
        g=self.gamma(etap)
        c=self.com(etap)
        value=-(n0*k/k0)*2.0*g**2/(1+g**2)
        return value

    def sc2(self,etap):
        n,n0=self.n,self.n0
        k,k0,kp=self.k,self.k0,self.kp
        g=self.gamma(etap)
        c=self.com(etap)
        value=(n0*k/k0)*2.0*g/(1+g**2)
        return value

    def sc3(self,etap):
        n,n0=self.n,self.n0
        k,kp=self.k,self.kp
        d=self.delta(etap)
        c=self.com(etap)
        value=2*(n/n0)*(k/kp)*d*c/(1+d**2)
        return value

    def sc4(self,etap):
        n,n0=self.n,self.n0
        k,kp=self.k,self.kp
        d=self.delta(etap)
        c=self.com(etap)
        value=2*(n/n0)*(k/kp)*c*d**2/(1+d**2)
        return value

    def sc5(self,etap):
        n,n0=self.n,self.n0
        k,k0,kp=self.k,self.k0,self.kp
        eta=self.Eta(etap)
        q=k*sin(eta)
        e=self.epsilon(etap)
        c=self.com(etap)
        value=(n/n0)*(q/kp)*(k/k0)*2/(1+e**2)
        return value

    def sc6(self,etap):
        n,n0=self.n,self.n0
        k,k0,kp=self.k,self.k0,self.kp
        eta=self.Eta(etap)
        q=k*sin(eta)
        e=self.epsilon(etap)
        c=self.com(etap)
        value=-(n/n0)*(q/kp)*(k/k0)*2*e/(1+e**2)
        return value

    # The electric field functions
    
    def Epperp(self,etap):
        n=self.n
        n0=self.n0
        np=self.np
        k0=self.k0
        k=self.k
        kp=self.kp

        eta=self.Eta(etap)
        eta0=self.Eta0(etap)
        q=k*sin(eta)
                        
        epperp=2.0*(q/kp)*(n/n0)*(k/k0)*sin(eta)*cos(eta)\
                /sin(eta+eta0)/cos(eta0-eta)
        return epperp

    def Eppar(self,etap):
        n=self.n
        n0=self.n0
        k=self.k
        k0=self.k0
        kp=self.kp

        eta=self.Eta(etap)
        eta0=self.Eta0(etap)

        eppar=2.0*(n/n0)*(k/kp)*cos(eta)*cos(eta0)*sin(eta)\
               /sin(eta+eta0)/cos(eta0-eta)
        return eppar
        
    def Espar(self,etap):
        n=self.n
        n0=self.n0
        k=self.k
        k0=self.k0
        kp=self.kp

        eta=self.Eta(etap)
        eta0=self.Eta0(etap)
    
        espar=-2.0*n*(k/k0)*cos(eta)*sin(eta)/sin(eta+eta0)
        return espar

    # Bessel functions of appropriate arguments
    
    def J0(self,etap):
        kp=self.kp
        rho=self.GetRho()
        j0=jn(0,kp*rho*etap)
        return j0

    def J1(self,etap):
        kp=self.kp
        rho=self.GetRho()
        j1=jn(1,kp*rho*etap)
        return j1
    
    def J2(self,etap):
        kp=self.kp
        rho=self.GetRho()
        j2=jn(2,kp*rho*etap)
        return j2
    
    # Calculation of the intensity
    
    def Val0_subreal(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))*self.J0(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)            
        return value

    def Val0_subimag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))*self.J0(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)            
        return value

    def Val0_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))*\
                       self.J0(etap)*cos(self.kp*self.deltaz*cos(eta))\
                       -etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))*\
                       self.J0(etap)*sin(self.kp*self.deltaz*cos(eta))
            
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val0_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))*\
                       self.J0(etap)*cos(self.kp*self.deltaz*cos(eta))\
                       +etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))*\
                       self.J0(etap)*sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val1_subreal(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.Epperp(etap)*self.J1(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)            
        return value

    def Val1_subimag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.Epperp(etap)*self.J1(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)            
        return value

    def Val1_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.sc5(etap)*self.J1(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       -etap/sqrt(cos(eta))*self.sc6(etap)*self.J1(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val1_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.sc6(etap)*self.J1(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       +etap/sqrt(cos(eta))*self.sc5(etap)*self.J1(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val2_subreal(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))*self.J2(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)            
        return value

    def Val2_subimag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))*self.J2(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)            
        return value

    def Val2_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))*self.J2(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       -etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))*self.J2(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val2_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))*self.J2(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       +etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))*self.J2(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    # Functions for the cumulant expansion
    
    def Bessel0Cum_subreal(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b0subreal_norm[self.focusindex]
        mom1=self.b0subreal_mean[self.focusindex]
        mom2=self.b0subreal_var[self.focusindex]

        arg=self.kp*rho*mom1
        j0=jn(0,arg)
        j1=jn(1,arg)
        j2=jn(2,arg)
        value=norm*(j0+0.5*(rho*self.kp)**2*mom2*(j2-j1/arg))
        return value

    def Bessel0Cum_subimag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b0subimag_norm[self.focusindex]
        mom1=self.b0subimag_mean[self.focusindex]
        mom2=self.b0subimag_var[self.focusindex]

        if mom1!=0.0:
            arg=self.kp*rho*mom1
            j0=jn(0,arg)
            j1=jn(1,arg)
            j2=jn(2,arg)
            value=norm*(j0+0.5*(rho*self.kp)**2*mom2*(j2-j1/arg))
        else:
            value=0.0
        return value

    def Bessel0Cum_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b0real_norm[self.focusindex]
        mom1=self.b0real_mean[self.focusindex]
        mom2=self.b0real_var[self.focusindex]

        arg=self.kp*rho*mom1
        value=norm*(jn(0,arg)+0.5*(rho*self.kp)**2*mom2*(jn(2,arg)-jn(1,arg)/arg))
        return value

    def Bessel0Cum_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b0imag_norm[self.focusindex]
        mom1=self.b0imag_mean[self.focusindex]
        mom2=self.b0imag_var[self.focusindex]

        arg=self.kp*rho*mom1
        value=norm*(jn(0,arg)+0.5*(rho*self.kp)**2*mom2*(jn(2,arg)-jn(1,arg)/arg))
        return value
    
    def Bessel1Cum_subreal(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b1subreal_norm[self.focusindex]
        mom1=self.b1subreal_mean[self.focusindex]
        mom2=self.b1subreal_var[self.focusindex]

        arg=self.kp*rho*mom1
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel1Cum_subimag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b1subimag_norm[self.focusindex]
        mom1=self.b1subimag_mean[self.focusindex]
        mom2=self.b1subimag_var[self.focusindex]

        if mom1!=0.0:
            arg=self.kp*rho*mom1
            j1=jn(1,arg)
            value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        else:
            value=0.0
        return value

    def Bessel1Cum_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b1real_norm[self.focusindex]
        mom1=self.b1real_mean[self.focusindex]
        mom2=self.b1real_var[self.focusindex]

        arg=self.kp*rho*mom1
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel1Cum_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b1imag_norm[self.focusindex]
        mom1=self.b1imag_mean[self.focusindex]
        mom2=self.b1imag_var[self.focusindex]

        arg=self.kp*rho*mom1
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel2Cum_subreal(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b2subreal_norm[self.focusindex]
        mom1=self.b2subreal_mean[self.focusindex]
        mom2=self.b2subreal_var[self.focusindex]

        arg=self.kp*rho*mom1
        j2=jn(2,arg)
        value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        return value

    def Bessel2Cum_subimag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b2subimag_norm[self.focusindex]
        mom1=self.b2subimag_mean[self.focusindex]
        mom2=self.b2subimag_var[self.focusindex]

        if mom1!=0.0:
            arg=self.kp*rho*mom1
            j2=jn(2,arg)
            value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                        (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        else:
            value=0.0
        return value

    def Bessel2Cum_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        kp=self.kp
        rho=self.GetRho()

        norm=self.b2real_norm[self.focusindex]
        mom1=self.b2real_mean[self.focusindex]
        mom2=self.b2real_var[self.focusindex]

        arg=self.kp*rho*mom1
        j2=jn(2,arg)
        value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        return value

    def Bessel2Cum_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed
        M=self.M
        np=self.np
        n=self.n

        rho=self.GetRho()

        norm=self.b2imag_norm[self.focusindex]
        mom1=self.b2imag_mean[self.focusindex]
        mom2=self.b2imag_var[self.focusindex]

        arg=self.kp*rho*mom1
        j2=jn(2,arg)
        value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        return value

    # Functions for calculation of the normalization constants
    
    def Intensity_norm_par(self,rho):
        self.SetRho(rho)
        par= ((self.Val0_subreal()+self.Val0_real())**2+(self.Val0_subimag()+self.Val0_imag())**2\
              +(self.Val2_subreal()+self.Val2_real())**2+(self.Val2_subimag()+self.Val2_imag())**2)\
              /4.0
        return par

    def Intensity_norm_vert(self,rho):
        self.SetRho(rho)
        vert=((self.Val1_subreal()+self.Val1_real())**2+(self.Val1_subimag()+self.Val1_imag())**2)
        return vert

    def NormIntegrand_par(self,r):
        rho=self.M*r
        value=zeros(len(rho))
        for i in range(len(rho)):
            value[i]=r[i]*self.Intensity_norm_par(rho[i])
        value*=(2*pi*self.M**2)
        return value

    def NormIntegrand_vert(self,r):
        rho=self.M*r
        value=zeros(len(rho))
        for i in range(len(rho)):
            value[i]=r[i]*self.Intensity_norm_vert(rho[i])
        value*=(2*pi*self.M**2)
        return value

    def Normalization(self):
        print("Calculating normalization constants...")
        norm_par=qromb(self.NormIntegrand_par,0.0,20000,1e-4)
        norm_vert=qromb(self.NormIntegrand_vert,0.0,20000,1e-4)             
        return (norm_par,norm_vert)
    
    # Calculation of the PSF

    def Intensity_exact(self,rho,phip,alpha,beta,deltaz):

        self.focusindex=find_nearest(self.focusvals,deltaz)
        self.deltaz=deltaz
        
        self.SetRho(rho)

        par= sin(beta)**2/4.0*\
             ((self.Val0_subreal()+self.Val0_real())**2+(self.Val0_subimag()+self.Val0_imag())**2\
              +(self.Val2_subreal()+self.Val2_real())**2+(self.Val2_subimag()+self.Val2_imag())**2\
               -2*cos(2*(phip-alpha))*((self.Val0_subreal()+self.Val0_real())*\
              (self.Val2_subreal()+self.Val2_real())\
               +(self.Val0_subimag()+self.Val0_imag())*(self.Val2_subimag()+self.Val2_imag())))

        mix=sin(beta)*cos(beta)*cos(phip-alpha)*\
             ((self.Val1_subreal()+self.Val1_real())*(self.Val2_subimag()+self.Val2_imag())-
              (self.Val2_subreal()+self.Val2_real())*(self.Val1_subimag()+self.Val1_imag())-
              (self.Val1_subreal()+self.Val1_real())*(self.Val0_subimag()+self.Val0_imag())+
              (self.Val0_subreal()+self.Val0_real())*(self.Val1_subimag()+self.Val1_imag()))

        vert=cos(beta)**2*\
              ((self.Val1_subreal()+self.Val1_real())**2+(self.Val1_subimag()+self.Val1_imag())**2)
        
        value=par+mix+vert

        # Normalization
        value/=(sin(beta)**2*self.norm[0]+cos(beta)**2*self.norm[1])
        return value


    def Intensity_approx(self,rho,phip,alpha,beta,deltaz):

        self.focusindex=find_nearest(self.focusvals,deltaz)
        self.deltaz=deltaz
        
        self.SetRho(rho)

        par= sin(beta)**2/4.0*\
             ((self.Bessel0Cum_subreal()+self.Bessel0Cum_real())**2+\
              (self.Bessel0Cum_subimag()+ self.Bessel0Cum_imag())**2\
              +(self.Bessel2Cum_subreal()+self.Bessel2Cum_real())**2+\
              (self.Bessel2Cum_subimag()+self.Bessel2Cum_imag())**2\
               -2*cos(2*(phip-alpha))*((self.Bessel0Cum_subreal()+self.Bessel0Cum_real())*\
              (self.Bessel2Cum_subreal()+self.Bessel2Cum_real())\
               +(self.Bessel0Cum_subimag()+self.Bessel0Cum_imag())*\
                                       (self.Bessel2Cum_subimag()+self.Bessel2Cum_imag())))

        mix=sin(beta)*cos(beta)*cos(phip-alpha)*\
             ((self.Bessel1Cum_subreal()+self.Bessel1Cum_real())*\
              (self.Bessel2Cum_subimag()+self.Bessel2Cum_imag())-\
              (self.Bessel2Cum_subreal()+self.Bessel2Cum_real())*\
              (self.Bessel1Cum_subimag()+self.Bessel1Cum_imag())-\
              (self.Bessel1Cum_subreal()+self.Bessel1Cum_real())*\
              (self.Bessel0Cum_subimag()+self.Bessel0Cum_imag())+\
              (self.Bessel0Cum_subreal()+self.Bessel0Cum_real())*\
              (self.Bessel1Cum_subimag()+self.Bessel1Cum_imag()))

        vert=cos(beta)**2*\
              ((self.Bessel1Cum_subreal()+self.Bessel1Cum_real())**2+\
               (self.Bessel1Cum_subimag()+self.Bessel1Cum_imag())**2)
        
        value=par+mix+vert

        # Normalization
        value/=(sin(beta)**2*self.norm[0]+cos(beta)**2*self.norm[1])
        return value

    def PSF_approx(self,x,y,alpha,beta,deltaz):
        # Transform alpha to define optical axis in direction from camera towards objective
        alpha+=pi
        self.deltaz=deltaz
        r=sqrt(x**2+y**2)
        if (x<0.0): phip=pi-arctan(-y/x)
        else: phip=arctan(y/x)
        rho=self.M*r
        value=self.M**2*self.Intensity_approx(rho,phip,alpha,beta,deltaz)
        return value

    def PSF_exact(self,x,y,alpha,beta,deltaz):
        # Transform alpha to define optical axis in direction from camera towards objective
        alpha+=pi
        self.deltaz=deltaz
        r=sqrt(x**2+y**2)
        if (x<0.0): phip=pi-arctan(-y/x)
        else: phip=arctan(y/x)
        rho=self.M*r
        value=self.M**2*self.Intensity_exact(rho,phip,alpha,beta,deltaz)
        return value

class LogLikelihood:
    """ Class defining the log-likelihood function maximized in MLE."""

    def __init__(self, counts, psf_generator):

        self.counts = counts
        self.psf_generator = psf_generator
        
        #self.pinit = pinit

    def Value(self,x):
        
        counts=self.counts
        # Transform the reparametrised theta and phi in the original theta and phi here
        #phi_conv = np.arctan2(x[1], x[0])%(2*np.pi)  # Reconvert x[0], x[1] to phi check the order of 0 and 1 needs to be the same as est
        #theta_conv = 0.5 * np.arccos(x[2])%(np.pi/2)  # Reconvert x[2] to theta

        # Alternative reparametrisation conversion
        phi_conv = np.arctan2(x[1], x[0])%(2*np.pi)
        theta_conv = np.arccos(np.clip(x[2], -1, 1))
        
        mux_nm = x[3]
        muy_nm = x[4]
        n_photons = x[5]

        # Calculate log-likelihood
        model_image = self.psf_generator(phi_conv, theta_conv, mux_nm, muy_nm, n_photons) # add verbose=True to activate it
        
        # Check for zero or negative values
        if np.any(n_photons <= 0):
            raise ValueError("n_photons must be positive! Received: {}".format(n_photons))
            #print("Warning: n_photons in model_image contains negative values!")
            #return -np.inf  
         
        ln_factorial = gammaln(counts + 1)
        log_likelihood = np.sum(counts * np.log(model_image) - model_image - ln_factorial)

        return -log_likelihood

class MLEwT:
    """
    Estimates the center coordinates (x and y) and the orientation (theta and phi)
    of a fixed dipole using MLEwT.

    Input:
    -------------------------------------------------------------------------

    wl (float)      : Wavelength [nm] of emission peak in buffer
    a (float)       : Width of pixels [nm] (assumed small)
    M (float)       : Magnification of the objective
    NA (float)      : Numerical aperture of the objective
    n (float)       : Refractive index of glass/immersion oil
    n0 (float)      : Refractive index of buffer
    alpha (float)   : Inverse gain of the EMCCD chip
    Sfloor (float)  : Constant offset of the EMCCD output
    sigma (float)   : Width of the noise distribution in the EMCCD output
    initvals (array): Array of length 6 of initial values for phi, theta, mux_nm, muy_nm,n_photons
    initpix (array) : Array of length 2 of initial values for the center pixel (ypixel,xpixel)
    deltapix (int)  : The half width of the array to be analyzed

    Functions:
    -------------------------------------------------------------------------
    Estimate (array) : Takes a full pixel array and uses MLEwT to return an
    array of estimates for x,y,b,N,theta,phi where b is the number of background
    photons per pixel and N is the photon number.

    Kim I. Mortensen
    """

    def __init__(self,wl,a,M,NA,n,np,n0,initvals,initpix,deltapix,Sfloor,alpha,sigma):

        # Store user input
        self.wl=wl
        self.a=a
        self.M=M
        self.NA=NA
        self.n=n
        self.np=np
        self.n0=n0

        self.Sfloor=Sfloor
        self.alpha=alpha
        self.sigma=sigma

        self.deltapix=deltapix

        self.initvals=initvals
        self.initpix=initpix

    def set_verbose(self, value):
        self._verbose = bool(value)

    def set_observer(self, value):
        assert hasattr(value, 'append')
        self._observer = value

    # Define the gradient (numerical derivative)
    def gradient(self, params):
        grad = np.zeros(len(params))  # Ensure it's a 1D array with length matching the number of params
        epsilon = 1e-6  # Perturbation value for finite differences
        for i in range(len(params)):
            perturbed = params.copy()
            perturbed[i] += epsilon
            grad[i] = (self.ll.Value(perturbed) - self.ll.Value(params)) / epsilon
        return grad

    def hessian(self, params):
        n_params = len(params)
        hess = np.zeros((n_params, n_params))  # Should be a 2D square matrix
        epsilon = 1e-6
        for i in range(n_params):
            for j in range(n_params):
                perturbed1 = params.copy()
                perturbed2 = params.copy()
                perturbed1[i] += epsilon
                perturbed2[j] += epsilon
                hess[i, j] = (self.ll.Value(perturbed1) + self.ll.Value(perturbed2) - 2 * self.ll.Value(params)) / (epsilon ** 2)
        return hess

    def Estimate(self,datamatrix):
        
        # Entire data matrix
        counts=datamatrix
        
        # Transformation of initial values
        if not isinstance(self.initvals, (list, np.ndarray)):
            raise ValueError(f"Expected initvals to be a list or np.ndarray of 6 elements, but got {self.initvals} {type(self.initvals)}")
        
        # Reparametrise phi -> sin(phi), cos(phi) this adds a parameter to the fitting
                       #theta -> cos(2*theta)
        pinit=zeros(6)
        pinit[0] = np.cos(self.initvals[0] % (2 * np.pi))       #self.initvals[0] phi reparametrised in cos(phi)
        pinit[1] = np.sin(self.initvals[0] % (2 * np.pi))       #self.initvals[1] was theta, now reparametrised phi in sin(phi)
        pinit[2] = np.cos(2 * self.initvals[1] % (np.pi / 2))   #reparametrised theta in cos(2theta)
        pinit[3] = self.initvals[2]               # mux_nm
        pinit[4] = self.initvals[3]               # muy_nm
        pinit[5] = self.initvals[4]               #sqrt(self.initvals[4]) # n_photons

        #print(f"pinit is: {pinit}")
        
        
        # Create instance of LogLikelihood object
        self.ll=LogLikelihood(counts, self.psf_generator)
        best_result = None
        best_ll_value = np.inf

        # Perform maximization of the log-likelihood using Powell's method
        #xopt, fopt, xi, direc, iter, funcalls, warnflag=\
        bounds = [(-1, 1), (-1, 1), (-1, 1), (-433, 433), (-433, 433), (1e-6, None)]

        options = {
             'maxiter': 100000000,  # Maximum number of iterations
             'disp': True,          # Enable display of progress
        }
        pinit[3] = np.clip(pinit[2], -433, 433)  # Ensure x (mux_nm) is within bounds
        pinit[4] = np.clip(pinit[3], -433, 433)  # Ensure y (muy_nm) is within bounds
        self.initvals[2] = np.clip(self.initvals[2], -433, 433)
        self.initvals[3] = np.clip(self.initvals[3], -433, 433)

        #res=fmin_powell(ll.Value,pinit, ftol=0.0001,maxiter=15,full_output=1, disp=False)
        res1 = minimize(self.ll.Value, pinit, method='Powell', bounds=bounds, options=options)
        first_attempt = res1.x
        first_cost = res1.fun

        # Check if theta is close to 90 degrees
        condition_value = abs(0.5 * np.arccos(first_attempt[2]) - np.pi / 2)
        print("condition_value is: ", condition_value)

        condition1 = condition_value >= 0.1 * np.pi / 180
        print("condition1 is: ", condition1)
        #condition1 = abs(0.5 * np.arccos(first_attempt[2]) - np.pi / 2) >= 0.1 * np.pi / 180

        if condition1:  # Good result, no need for additional attempts
            best_result = first_attempt
        else:
            # Prepare attempts
            attempts = [first_attempt]
            costs = [first_cost]

            # Second attempt (flip phi)
            second_start = pinit.copy()
            second_start[0] *= -1
            second_start[1] *= -1
            res2 = minimize(self.ll.Value, second_start, method='Powell', bounds=bounds, options=options)
            attempts.append(res2.x)
            costs.append(res2.fun)

            # Third attempt (flip theta)
            third_start = pinit.copy()
            third_start[2] *= -1
            res3 = minimize(self.ll.Value, third_start, method='Powell', bounds=bounds, options=options)
            attempts.append(res3.x)
            costs.append(res3.fun)

            # Fourth attempt (flip both theta and phi)
            fourth_start = pinit.copy()
            fourth_start[0] *= -1
            fourth_start[1] *= -1
            fourth_start[2] *= -1
            res4 = minimize(self.ll.Value, fourth_start, method='Powell', bounds=bounds, options=options)
            attempts.append(res4.x)
            costs.append(res4.fun)

            # Select the best attempt
            best_index = np.argmin(costs)
            best_result = attempts[best_index]

        # Debugging gradient and hessian
        grad = self.gradient(best_result)
        hess = self.hessian(best_result)

        print(f"Gradient: {grad}, Shape: {grad.shape}")
        print(f"Hessian: {hess}, Shape: {hess.shape}")

        try:
            res_refined = minimize(self.ll.Value, best_result, method='trust-constr', jac=self.gradient(best_result), hess=self.hessian(best_result), bounds=bounds, options=options)
            best_result = res_refined.x
        except Exception as e:
            print(f"Error during refinement: {e}")

        # Refinement using trust-constr if needed
        #options = {'disp': self._verbose}
        #res_refined = minimize(self.ll.Value, best_result, method='trust-constr', jac=self.gradient(best_result), hess=self.hessian(best_result), bounds=bounds, options=options)
        #best_result = res_refined.x

        # Transform results back to original parameters
        est = best_result
        self.phi_conv = np.arctan2(est[1], est[0])
        self.theta_conv = np.arccos(est[2])
        self.mux_nm = est[3]
        self.muy_nm = est[4]
        n_photons = est[5]

        results = [self.phi_conv, self.theta_conv, self.mux_nm, self.muy_nm, n_photons]
        return results
        
        # Convert estimates
        #est[2]=est[2]**2
        #est[3]=est[3]**2
        #print(f"est is: {est}")

        # Calculate covariance matrix of estimates of position coordinates and angles
        #covar=MLEwTcovar(self.a,self.deltapix*2,self.wl,self.n,
        #                 self.n0,self.M,self.NA, self._norm_file)

        #covarmatrix=covar.CovarianceMatrix(est[3],est[2],\
        #                             array([self.mux,self.muy]),est[4],est[5])

        # Add escess noise
        #covarmatrix*=2.0

        #errorbars=sqrt(diag(covarmatrix))

        #self._observer.append([
        #    self.a*xpix+est[0], errorbars[0],
        #    self.a*ypix+est[1], errorbars[1],
        #    est[2], est[3],
        #    est[4], errorbars[2],
        #    est[5], errorbars[3],
        #    ])

        #if self._verbose:
        #    print("\nx coordinate [nm] = ", around(self.a*xpix+est[0],1),'+/-',around(errorbars[0],1))
        #    print("y coordinate [nm] = ", around(self.a*ypix+est[1],1),'+/-',around(errorbars[1],1))
        #    print("azimuthal angle [rad] = ", around(est[4],2),'+/-',around(errorbars[2],3))
        #    print("polar angle [rad] = ", around(est[5],2),'+/-',around(errorbars[3],3))
        #    print('')
        #    print(covarmatrix)

        #return est

class MLEwTcovar:
    """
    Calculates the covariance matrix for the estimated parameters in MLEwT.

    The PSF is the distribution of photons in the image plane from
    a fixed dipole emitter close to the coverslip surface,
    when imaged by a large aperture objective.

    Input:
    -------------------------------------------------------------------------

    a (float)  : Width of pixels [nm] (assumed small)
    npix (int) : Number of analyzed pixels along one dimension
    wl (float) : Wavelength [nm] of emission peak in buffer.
    NA (float) : Numerical aperture of the objective
    n (float) : Refractive index of glass/immersion oil
    n0 (float) : Refractive index of buffer (water)
    M (float) : Magnification of the objective

    Functions:
    -------------------------------------------------------------------------
    CovarianceMatrix (float) : Takes the photon number, number of background photons
    per pixel, the center coordinates (array([x,y])), the azimuthal angle (theta), and the
    polar angle (phi) of the probe as arguments. Returns the covariance matrix for the estimated
    parameters x,y,theta, and phi.

    Kim I. Mortensen
    """

    def __init__(self,a,npix,wl,n,n0,M,NA, norm_file):

        self.a=a
        self.npix=npix
        self.wl=wl
        self.NA=NA
        self.n=n
        self.n0=n0
        self.M=M


        # Define instance of dipole PSF
        self.dip=dipdistr(self.wl,self.n,self.n0,self.M,self.NA, norm_file)

    #def Probability(self,mu,theta,phi):

    #    a=self.a
    #    npix=self.npix

        # Center pixel coordinates in one dimension
    #    posvec=arange(-(npix-1.0)/2.0,npix/2.0,1.0)*a
        
        # Matrix of distances from PSF center to pixel centers
        #distmat=zeros((npix,npix))
        # Matrix of expected photon counts
    #    p=zeros((npix,npix))

        # -- Calculate expected photon values in small pixels
        #    These are calculated by approximating the integral of
        #    the PSF over a pixel to leading order in the pixel width.

    #    for i in range(npix):
    #        for j in range(npix):
    #            x=posvec[i]-mu[0]
    #            y=posvec[j]-mu[1]

    #            p[j,i]=self.dip.PSF_approx(phi,theta,x,y)

        # Expected photon counts (using 1. order approximation over small pixels)
    #    p*=(a**2)
    #    return p

    def Derivatives(self,mu,theta,phi):

        delta=1e-6
        f0=(self.Probability(mu+array([delta,0.0]),theta,phi)\
             -self.Probability(mu-array([delta,0.0]),theta,phi))/(2*delta)
        f1=(self.Probability(mu+array([0.0,delta]),theta,phi)\
             -self.Probability(mu-array([0.0,delta]),theta,phi))/(2*delta)
        delta=1e-8
        f2=(self.Probability(mu,theta+delta,phi)\
             -self.Probability(mu,theta-delta,phi))/(2*delta)
        f3=(self.Probability(mu,theta,phi+delta)\
             -self.Probability(mu,theta,phi-delta))/(2*delta)

        return (f0,f1,f2,f3)

    def FisherMatrix(self,N,b,mu,theta,phi):

        p=self.Probability(mu,theta,phi)
        f1,f2,f3,f4=self.Derivatives(mu,theta,phi)

        I=zeros((4,4))

        denom=p+1.0*b/N

        I[0,0]=sum(ravel(f1**2/denom))
        I[0,1]=I[1,0]=sum(ravel(f1*f2/denom))
        I[0,2]=I[2,0]=sum(ravel(f1*f3/denom))
        I[0,3]=I[3,0]=sum(ravel(f1*f4/denom))
        I[1,1]=sum(ravel(f2**2/denom))
        I[1,2]=I[2,1]=sum(ravel(f2*f3/denom))
        I[1,3]=I[3,1]=sum(ravel(f2*f4/denom))
        I[2,2]=sum(ravel(f3**2/denom))
        I[2,3]=I[3,2]=sum(ravel(f3*f4/denom))
        I[3,3]=sum(ravel(f4**2/denom))

        I*=N
        return I

    def CovarianceMatrix(self,N,b,mu,theta,phi):
        return inv(self.FisherMatrix(N,b,mu,theta,phi))
