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
# dave feb 2025 - particle swarm
from pyswarm import pso

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
    
    """ Class defining the log-likelihood function maximized in MLE.

        Class initialization

        Input:
        ------------------------------------------------------------
        pw (float)      : Pixel width
        wl (float)      : Peak emission wavelength
        NA (float)      : Numerical aperture of objective
        n (float)       : Refractive index of immersion oil/optics
        n0 (float)      : Refractive index of buffer
        M (float)       : Magnification of composite microscope
        pinit (array)   : Array (len=6) of initial parameter values (mux,muy,n,N,theta,phi)
        Sfloor (float)  : Constant offset of the EMCCD camera
        alpha (float)   : Inverse gain of the EMCCD camera
        sigma (float)   : Width of the readout-noise distribution of the EMCCD camera

        Functions:
        ------------------------------------------------------------
        Value           : Calculates the negative log-likelihood
    """

    def __init__(self,counts,pw,wl,NA,n,n0,M,pinit,Sfloor,alpha,sigma):

        self.counts=counts
        self.pw=pw
        self.wl=wl
        self.NA=NA
        self.n=n
        self.n0=n0
        self.M=M

        self.Sfloor=Sfloor
        self.alpha=alpha
        self.sigma=sigma
        
        self.pinit=pinit
        
        self.npix=shape(counts)[0]

        # Pixel coordinate vector
        self.posvec=arange(-(self.npix-1.0)/2.0,(self.npix)/2.0,1.0)*pw

        # Create instance of dipole PSF
        self.dd=dipdistr(self.wl,self.n,self.n0,self.M,self.NA)
    
    # def Value(self,x):
    #
    #     """ Calculation of the negative log-likelihood
    #
    #         Input:
    #         ------------------------------------------------------------
    #         x (array)       : Array of current parameter values
    #
    #         Output:
    #         ------------------------------------------------------------
    #         value (float)   : The negative log-likelihood
    #     """
    #
    #     counts=self.counts
    #     npix=self.npix
    #     posvec=self.posvec
    #
    #     Sfloor=self.Sfloor
    #     alpha=self.alpha
    #     sigma=self.sigma
    #
    #     pij=zeros((npix,npix))
    #
    #     params=zeros(7)
    #
    #     # Assign parameters their current values
    #     mux,muy,b,N,phi,theta,dz=x
    #
    #     # Convert parameters
    #     b=b**2
    #     N=N**2
    #
    #     # Calculate probabilities for all pixels
    #     for i in range(int(npix)):
    #         for j in range(int(npix)):
    #             pij[j,i]=self.dd.PSF_approx(posvec[i]-mux,posvec[j]-muy,phi,theta,dz)
    #     pij*=(self.pw**2)
    #
    #     # Subtract noise floor
    #     effcounts=counts-Sfloor
    #
    #     # Calculate log-likelihood
    #     value=0.0
    #     for i in range(int(npix)):
    #         for j in range(int(npix)):
    #             eta=N*pij[j,i]+b
    #
    #             f0=alpha*exp(-eta)*eta
    #             fp0=f0*0.5*alpha*(eta-2)
    #
    #             cij=effcounts[j,i]
    #
    #             # Functions for approximate convolution of the EMCCD readout-noise distribution
    #             # with the signal distribution from the amplification
    #             conv0=0.5*(1+erf(cij/(sqrt(2*sigma**2))))
    #             conv1=sigma*exp(-cij**2/(2*sigma**2))/sqrt(2*pi)+cij*conv0
    #             temp=(f0*conv0+fp0*conv1+exp(-eta)*gauss(cij,sigma))
    #
    #             if (cij>0.0):
    #                 nij=alpha*cij
    #                 # Use log-transform for large arguments
    #                 if eta*nij>10**5:
    #                     transform=0.5*log(alpha*eta/cij)-nij-eta+2*sqrt(eta*nij)\
    #                                -log(2*sqrt(pi)*(eta*nij)**0.25)
    #                     temp+=(exp(transform)-f0-fp0*cij)
    #                 else:
    #                     temp+=(sqrt(alpha*eta/cij)*exp(-nij-eta)*i1(2*sqrt(eta*nij))\
    #                         -f0-fp0*cij)
    #
    #             # dave feb 2025 - temp sometimes hits 0 due to rounding errors
    #             if temp <= 0:
    #                 temp = 1e-6
    #
    #             value+=log(temp)
    #
    #     # The negative log-likelihood is to be minimized
    #     value*=-1.0
    #     # dave nov 2024 - commented out this below
    #     # print("{:7.5f} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f}".format(value, mux, muy, b, N, theta, phi, dz))
    #     return value

    # dave feb 2025 - fix photon number, background noise, and z
    def Value(self, x):

        """ Calculation of the negative log-likelihood

            Input:
            ------------------------------------------------------------
            x (array)       : Array of current parameter values

            Output:
            ------------------------------------------------------------
            value (float)   : The negative log-likelihood
        """

        counts = self.counts
        npix = self.npix
        posvec = self.posvec

        Sfloor = self.Sfloor
        alpha = self.alpha
        sigma = self.sigma

        pij = zeros((npix, npix))

        params = zeros(8)

        # Assign parameters their current values
        mux, muy, b, N, cosphi, sinphi, costheta, dz = x

        # dave feb 2025 - manually put in bounds on costheta etc.
        cosphi = max(min(cosphi, 1), -1)
        sinphi = max(min(cosphi, 1), -1)
        costheta = max(min(cosphi, 1), -1)
        # dave feb 2025 - reparam
        phi = np.atan2(sinphi, cosphi)
        theta = 0.5*np.acos(costheta)


        # Convert parameters
        b = b ** 2
        N = N ** 2

        # Calculate probabilities for all pixels
        for i in range(int(npix)):
            for j in range(int(npix)):
                pij[j, i] = self.dd.PSF_approx(posvec[i] - mux, posvec[j] - muy, phi, theta, dz)
        pij *= (self.pw ** 2)

        # Subtract noise floor
        effcounts = counts - Sfloor

        # Calculate log-likelihood
        value = 0.0
        for i in range(int(npix)):
            for j in range(int(npix)):
                eta = N * pij[j, i] + b

                f0 = alpha * exp(-eta) * eta
                fp0 = f0 * 0.5 * alpha * (eta - 2)

                cij = effcounts[j, i]

                # Functions for approximate convolution of the EMCCD readout-noise distribution
                # with the signal distribution from the amplification
                conv0 = 0.5 * (1 + erf(cij / (sqrt(2 * sigma ** 2))))
                conv1 = sigma * exp(-cij ** 2 / (2 * sigma ** 2)) / sqrt(2 * pi) + cij * conv0
                temp = (f0 * conv0 + fp0 * conv1 + exp(-eta) * gauss(cij, sigma))

                if (cij > 0.0):
                    nij = alpha * cij
                    # Use log-transform for large arguments
                    if eta * nij > 10 ** 5:
                        transform = 0.5 * log(alpha * eta / cij) - nij - eta + 2 * sqrt(eta * nij) \
                                    - log(2 * sqrt(pi) * (eta * nij) ** 0.25)
                        temp += (exp(transform) - f0 - fp0 * cij)
                    else:
                        temp += (sqrt(alpha * eta / cij) * exp(-nij - eta) * i1(2 * sqrt(eta * nij)) \
                                 - f0 - fp0 * cij)

                # dave feb 2025 - temp sometimes hits 0 due to rounding errors
                if temp <= 0:
                    temp = 1e-6

                value += log(temp)

        # The negative log-likelihood is to be minimized
        value *= -1.0
        # dave nov 2024 - commented out this below
        # print("{:7.5f} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f}".format(value, mux, muy, b, N, theta, phi, dz))
        return value

    
class MLEwT:
    
    """ Maximum likelihood estimation of fluorophore center coordinates, orientation angles,
        and defocus.

        This class uses the function Estimate to perform maximum likelihood estimation of
        location and orientation of a single, fixed fluorophore and the objective's distance
        from the design focal plane using a theoretical dipole PSF as calculated in the module 
        dipdistr.py. The negative log-likelihood to be minimized is calculated using the 
        class LogLikelihood, while the covariance matrix is calculated using
        the class MLEwTcovar.py.

        Class initiation

        Input parameters:
        ------------------------------------------------------------ 
        wl (float)      : Peak emission wavelength
        a (float)       : Pixel width
        M (float)       : Magnification of composite microscope
        NA (float)      : Numerical aperture of objective
        n (float)       : Refractive index of immersion oil/optics
        np (float)      : Refractive index of imaging medium (typically that of air, i.e. np=1.0)
        n0 (float)      : Refractive index of buffer
        initvals (array): Array (len=7) of initial parameter values (mux,muy,b,N,phi,theta,dz)
        initpix (array) : Pixel coordinates of the center pixel in the area to be analyzed
        deltapix (int)  : Half-width of the pixel area to be analyzed
        Sfloor (float)  : Constant offset of the EMCCD camera
        alpha (float)   : Inverse gain of the EMCCD camera
        sigma (float)   : Width of the readout-noise distribution of the EMCCD camera

        Functions:
        ------------------------------------------------------------
        Estimate        : Performs maximum likelihood estimation on a spot for a single frame
        
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


    def Estimate(self,datamatrix):
        """ Performs maximum likelihood estimation from a diffraction-limited spot."""

        ypix=self.initpix[0]
        xpix=self.initpix[1]
        deltapix=self.deltapix

        # Entire data matrix
        counts=datamatrix

        # Extract pixel array around initial pixel
        counts=counts[int(ypix-deltapix):int(ypix+deltapix),int(xpix-deltapix):int(xpix+deltapix)]

        # dave feb 2025
        # reparameterise: theta -> cos(2*theta)
        # phi -> sin(phi), cos(phi)

        # Transformation of initial values
        # pinit=zeros(7)
        # pinit[0:2]=self.initvals[0:2] # x, y
        # pinit[2]=sqrt(self.initvals[2]) # b
        # pinit[3]=sqrt(self.initvals[3]) # N
        # pinit[4]=self.initvals[4] # phi
        # pinit[5]=self.initvals[5] # theta
        # pinit[6]=self.initvals[6] # dz
        pinit=zeros(8)
        pinit[0:2]=self.initvals[0:2] # x, y
        pinit[2]=sqrt(self.initvals[2]) # b
        pinit[3]=sqrt(self.initvals[3]) # N
        pinit[4]=np.cos(self.initvals[4] % 2*pi) # cos(phi)
        pinit[5]=np.sin(self.initvals[4] % 2*pi) # sin(phi)
        pinit[6]=np.cos(2*(self.initvals[5] % 2*pi)) # cos(theta)
        pinit[7]=self.initvals[6] # dz

        # Create instance of LogLikelihood object
        # ll=LogLikelihood(counts,self.a,self.wl,self.NA,self.n,self.n0,self.M,pinit,\
        #                  self.Sfloor,self.alpha,self.sigma)
        ll=LogLikelihood(counts,self.a,self.wl,self.NA,self.n,self.n0,self.M,pinit,\
                         self.Sfloor,self.alpha,self.sigma)

        # Perform maximization of the log-likelihood using Powell's method
        start_powell = time.time()
        # res=fmin_powell(ll.Value,pinit,ftol=0.000000001,maxiter=50,maxfun=500)

        bounds = [(-1000, 1000), (-1000, 1000), (None, None), (None, None), (-1, 1), (-1, 1), (-1, 1), (None, None)]

        options = {
            'maxiter': 1000,  # Maximum number of iterations
            # 'maxfev': 500,  # Maximum number of function evaluations
        }
        res=minimize(ll.Value,pinit,method='trust-constr',bounds=bounds, options=options)

        end_powell = time.time()
        elapsed_time_powell = end_powell - start_powell
        print(f"Time: {elapsed_time_powell:.4f} seconds to maximization log-likelihood")

        est=res.x#[0]
        # warnflag=res[5]

        # L=ll.Value(est)

        # Store position estimates relative to initial pixel
        self.mux=est[0]
        self.muy=est[1]

        # Convert estimates
        est[2]=est[2]**2
        est[3]=est[3]**2

        # Calculate covariance matrix
        asegment=self.a
        self.npix=2*self.deltapix
        covar=MLEwTcovar(self.a,self.npix,self.wl,self.NA,self.n,self.n0,self.M,asegment)

        # # dave jan 2025 - commenting out covar because of errors
        # covarmatrix=covar.CovarianceMatrix(est[3],est[2],est[0:2],est[4],est[5],est[6])
        #
        # # Account for excess noise
        # covarmatrix*=2.0
        #
        # errorbars=sqrt(diag(covarmatrix))
        #
        # print("\nx coordinate [nm] = ", around(est[0],3),'+/-',around(errorbars[0],1))
        # print("y coordinate [nm] = ", around(est[1],3),'+/-',around(errorbars[1],1))
        # print("azimuthal angle [rad] = ", around(est[4],2),'+/-',around(errorbars[2],3))
        # print("polar angle [rad] = ", around(est[5],2),'+/-',around(errorbars[3],3))

        # # dave nov 2024 - i also commented out the print statements above
        # print("covariance matrix = ")
        # print(covarmatrix)

        x_est = est[0]
        y_est = est[1]
        # undo the reparameterisation
        phi_est = atan2(est[5], est[4])# % (2 * np.pi)
        theta_est = 0.5*acos(est[6])# % (2 * np.pi)
        # cov_mat = covarmatrix

        print("\nx coordinate [nm] = ", around(x_est,3))
        print("y coordinate [nm] = ", around(y_est,3))
        print("azimuthal angle [rad] = ", around(phi_est,2))
        print("polar angle [rad] = ", around(theta_est,2))




        # # Define the radial and azimuthal range
        # r_vals = np.linspace(0, 29*51.2/2, 6)  # Radial values (in pixels or nm, adjust as needed)
        # azimuth_vals = np.linspace(0, 2 * np.pi, 6)  # Azimuth values (0 to 2*pi)
        #
        # # Initialize cost function values (Negative Log-Likelihood)
        # cost_matrix = np.zeros((len(r_vals), len(azimuth_vals)))
        #
        # # Compute the objective function (Negative Log-Likelihood) over radial and azimuthal coordinates
        # for i in range(len(r_vals)):
        #     for j in range(len(azimuth_vals)):
        #         # Convert (r, θ) to (x, y)
        #         x = r_vals[i] * np.cos(azimuth_vals[j])
        #         y = r_vals[i] * np.sin(azimuth_vals[j])
        #
        #
        #         # Transformation of initial values
        #         pinit=zeros(7)
        #         pinit[0:2]=i
        #         pinit[2]=sqrt(self.initvals[2])
        #         pinit[3]=sqrt(self.initvals[3])
        #         pinit[4]=j
        #         pinit[5]=self.initvals[5]
        #         pinit[6]=self.initvals[6]
        #
        #         # Create instance of LogLikelihood object
        #         ll=LogLikelihood(counts,self.a,self.wl,self.NA,self.n,self.n0,self.M,pinit,\
        #                          self.Sfloor,self.alpha,self.sigma)
        #
        #         # Perform maximization of the log-likelihood using Powell's method
        #         res=fmin_powell(ll.Value,pinit,ftol=0.0000001,maxiter=5000,full_output=1)#,maxfun=100000)
        #         est=res[0]
        #
        #         cost_matrix[i, j] = -ll.Value(est)
        #
        # # Plot the objective function as a 3D surface plot (or 2D contour plot)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        #
        # # Convert polar to Cartesian for surface plot
        # AZI, R = np.meshgrid(azimuth_vals, r_vals)
        # X, Y = R * np.cos(AZI), R * np.sin(AZI)
        #
        # # Create the 3D surface plot
        # ax.plot_surface(X, Y, cost_matrix, cmap='jet', edgecolor='none')
        # ax.set_xlabel('X Position')
        # ax.set_ylabel('Y Position')
        # ax.set_zlabel('Objective Function Value')
        # ax.set_title('Objective Function in Polar Coordinates (Inclination = 0)')
        #
        # plt.show()



        return x_est, y_est, theta_est, phi_est, 0#cov_mat
        # return


class MLEwTcovar:
    
    """ Class for calculating the covariance matrix of estimates obtained with MLEwT.

        The covariance matrix is calculated as the inverse Fisher information matrix.
        Probabilities and derivatives are calculated over a discrete array of imaginary, small
        pixels, to allow for a constant approximation of the PSF over such a pixel. To get the
        appropriate probabilities for the pixels used in the estimation, the probabilities of
        these small, imaginary pixels are then added.

        The calculation of this covariance matrix assumes that both center coordinates and
        both angles of orientation have been fitted.

        Input:
        ----------------------------------------------------------------------------
        a (float)       : Pixel width
        npix (int)      : Full width of the analyzed pixel array
        wl (float)      : Peak emission wavelength of the fluorophore
        NA (float)      : Numerical aperture of the objective
        n (float)       : Refractive index of immersion oil/optics
        n0 (float)      : Refractive index of the buffer (typically that of water)
        M (float)       : Magnification of the composite microscope
        asegment (float): The width of imaginary pixels in an array used for the calculation
                          of the covariance matrix. This should be chosen such that a constant
                          approximation of the PSF over such an imaginary pixel is adequate.

        Functions:
        ----------------------------------------------------------------------------
        Probability     : Calculates the array of pixel probabilities
        Derivatives     : Calculates array of numerical derivatives for the center coordinates
                          and orientation angles
        FisherMatrix    : Calculates Fisher's information matrix for given photon number
                          background level
        CovarianceMatrix: Calculates the covariance matrix of the estimate
        Variances       : Extracts the variances from the covariance matrix
                          
        """

    def __init__(self,a,npix,wl,NA,n,n0,M,asegment):

        self.a=a
        self.npix=npix
        self.wl=wl
        self.NA=NA
        self.n=n
        self.n0=n0
        self.M=M
        self.asegment=asegment

        # Define instance of dipole PSF
        self.dip=dipdistr(self.wl,self.n,self.n0,self.M,self.NA)

    def Probability(self,mu,phi,theta,deltaz):
        """ Calculates the array of probabilities over an array of pixels when the probe is
            located at mu=(x,y) and oriented with angles theta and phi."""

        a=self.a
        asegment=self.asegment
        npix=self.npix

        # Center pixel coordinates in one dimension
        posvec=arange(-(npix-1.0)/2.0,npix/2.0,1.0)*a
        # Matrix of distances from PSF center to pixel centers
        distmat=zeros((int(npix),int(npix)))
        # Matrix of expected photon counts 
        p=zeros((int(npix),int(npix)))

        # -- Calculate expected photon values in pixels
        #    These are calculated by approximating the integral of
        #    the PSF over a pixel to leading order in the pixel width.
        #    Large pixels are subdivided into smaller pixel segments.

        nsub=ceil(1.0*a/asegment)   # #-divisions

        if mod(nsub,2.0)==0.0:      # Make sure that #-divisions is uneven
            nsub+=1
            
        asub=a/nsub                 # Width of pixel segments

        if nsub>1: 
            for i in range(int(npix)):
                for j in range(int(npix)):
                    subrange=append(-1*arange((nsub+1)/2)[1:][::-1],arange((nsub+1)/2))
                    for k in subrange:
                        for l in subrange:
                            x=posvec[i]+l*asub-mu[0]
                            y=posvec[j]+k*asub-mu[1]
                            
                            # Expected photon counts (using 1. order approximation over small pixels)
                            p[j,i]+=(asub**2*self.dip.PSF_approx(x,y,phi,theta,deltaz))

        else: 
            for i in range(int(npix)):
                for j in range(int(npix)):
                    x=posvec[i]-mu[0]
                    y=posvec[j]-mu[1]
            
                    p[j,i]=self.dip.PSF_approx(x,y,phi,theta,deltaz)

            # Expected photon counts (using 1. order approximation over small pixels)
            p*=(a**2)

        return p

    def Derivatives(self,mu,phi,theta,deltaz):
        """ Calculates array of numerical derivatives for the center coordinates
            and orientation angles."""

        delta=1e-6
        f0=(self.Probability(mu+array([delta,0.0]),phi,theta,deltaz)\
             -self.Probability(mu-array([delta,0.0]),phi,theta,deltaz))/(2*delta)
        f1=(self.Probability(mu+array([0.0,delta]),phi,theta,deltaz)\
             -self.Probability(mu-array([0.0,delta]),phi,theta,deltaz))/(2*delta)
        delta=1e-8
        f2=(self.Probability(mu,phi+delta,theta,deltaz)\
             -self.Probability(mu,phi-delta,theta,deltaz))/(2*delta)
        f3=(self.Probability(mu,phi,theta+delta,deltaz)\
             -self.Probability(mu,phi,theta-delta,deltaz))/(2*delta)
        delta=2.0
        f4=(self.Probability(mu,phi,theta,deltaz+delta)\
             -self.Probability(mu,phi,theta,deltaz-delta))/(2*delta)

        return (f0,f1,f2,f3,f4)

    def FisherMatrix(self,N,b,mu,phi,theta,deltaz):
        """ Calculates Fisher's information matrix for given photon number N and background level
            b with the probe located at mu=(x,y) and oriented at angles theta and phi."""

        p=self.Probability(mu,phi,theta,deltaz)
        f1,f2,f3,f4,f5=self.Derivatives(mu,phi,theta,deltaz)

        I=zeros((5,5))

        denom=p+1.0*b/N

        I[0,0]=sum(ravel(f1**2/denom))
        I[0,1]=I[1,0]=sum(ravel(f1*f2/denom))
        I[0,2]=I[2,0]=sum(ravel(f1*f3/denom))
        I[0,3]=I[3,0]=sum(ravel(f1*f4/denom))
        I[0,4]=I[4,0]=sum(ravel(f1*f5/denom))
        I[1,1]=sum(ravel(f2**2/denom))
        I[1,2]=I[2,1]=sum(ravel(f2*f3/denom))
        I[1,3]=I[3,1]=sum(ravel(f2*f4/denom))
        I[1,4]=I[4,1]=sum(ravel(f2*f5/denom))
        I[2,2]=sum(ravel(f3**2/denom))
        I[2,3]=I[3,2]=sum(ravel(f3*f4/denom))
        I[2,4]=I[4,2]=sum(ravel(f3*f5/denom))
        I[3,3]=sum(ravel(f4**2/denom))
        I[3,4]=I[4,3]=sum(ravel(f4*f5/denom))
        I[4,4]=sum(ravel(f5**2/denom))

        I*=N
        return I

    def CovarianceMatrix(self,N,b,mu,phi,theta,deltaz):
        """ Calculates the covariance matrix of the estimate."""
        return inv(self.FisherMatrix(N,b,mu,phi,theta,deltaz))

    def Variances(self,N,b,mu,phi,theta,deltaz):
        """ Extracts the variances from the covariance matrix."""
        return diag(self.CovarianceMatrix(N,b,mu,phi,theta,deltaz))



if __name__=='__main__':

    close('all')

    # Experimental parameters
    wl=580.0
    n=1.52
    n0=1.33
    NA=1.49
    M=250.0
    
    deltax=32.0

    # PSF parameters
    N=10000.0
    b=5.0
    
    mu=0.1
    nu=0.1
    phi=2*pi/3.0
    theta=0.5

    # Distance from design focal plane
    deltaz=-30.0

    # EMCCD parameters
    inversegain=0.0089
    sigmanoise=12.6
    Sfloor=300.0
    gain=1.0/inversegain

    # Pixel array parameters
    npix=12
    deltapix=npix/2

    # Load data
    datamatrix=numpy.loadtxt('data.txt')

    # Initial guess
    initvals = array([mu,nu,b,N,phi,theta,deltaz])
    initpix=(deltapix,deltapix)

    # Create instance of MLEwT
    track=MLEwT(wl,deltax,M,NA,n,np,n0,initvals,initpix,deltapix,Sfloor,inversegain,sigmanoise)
    # Perform estimation
    track.Estimate(datamatrix)

    show()
