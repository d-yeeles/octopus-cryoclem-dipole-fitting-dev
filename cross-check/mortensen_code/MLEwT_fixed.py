# Python (python.org) code for the calculation of the theoretical PSF
# (Eq.(3) in Mortensen et al.) for fixed dipoles, maximum likelihood
# estimator for the position and orientation of the probe,
# and the calculation of the covariance matrix for the estimated
# parameters.
# Uses matplotlib, scipy.
#
#
# 16th of December 2009
#
# Kim I. Mortensen

# from scipy.optimize import *
# from pylab import *

import warnings
warnings.filterwarnings('error')

import numpy as np
import math
import os

from numpy import pi, shape, array, ravel, zeros, arange
from numpy import fabs, diag, around, exp, sqrt, log
from numpy import sin, cos
from numpy.linalg import inv
from scipy.special import i1, jn, gamma, erf, gammaln
from scipy.optimize import fmin_powell, minimize, BFGS
from scipy.optimize import SR1, BFGS
from numint import qromb
#from numint import qromb

def gauss(x,s):
    return exp(-x**2/(2*s**2))/sqrt(2*pi*s**2)

def fac(n):
    return gamma(n+1)

class dipdistr:

    """

    Calculates the theoretical point spread function (PSF) for fixed dipoles.

    The PSF is the distribution of photons in the image plane from
    a fixed dipole emitter close to the coverslip surface,
    when imaged by a large aperture objective.

    Input:
    -------------------------------------------------------------------------

    wl (float) : Wavelength [nm] of emission peak in buffer.
    NA (float) : Numerical aperture of the objective
    n (float) : Refractive index of glass/immersion oil
    n0 (float) : Refractive index of buffer (water)
    M (float) : Magnification of the objective

    Functions:
    -------------------------------------------------------------------------
    PSF_exact (float) : Takes the coordinates in the image and the angles of the
    probe. Returns the value of the exact PSF at that point.

    PSF_approx (float): Takes the coordinates in the image and the angles of the
    probe. Returns the value of the approximated PSF at that point.

    """

    def __init__(self, wavelen, n_objective, n_sample, magnification, NA, norm_file):

        self.wl = wavelen          #Wavelength in nm (emission peak in buffer)
        self.NA = NA          #Numerical aperture

        self.M = magnification            #Magnification
        self.n = n_objective            #Refractive index of glass/oil
        self.n_air = 1.0         #Refractive index of air in lab (vacuum)
        self.n0 = n_sample          #refractive index of sample medium (water)

        self.kp=2*np.pi/wavelen     #Wavevector amplitude in air (vacuum)
        self.k0=n_sample*self.kp  #Wavevector amplitude in buffer
        self.k=n_objective*self.kp    #Wavevector amplitude in glass

        #self.etapmed=n_sample/magnification   #Integration limits
        self.etapmed = 0.00609
        self.etapmax = 0.01009
        #self.etapmax=NA/magnification

        # Calculate or load normalization constants
        try:
            normdata=np.load(norm_file)

            try: shape(normdata)[1]
            except IndexError: normdata=array([normdata])

            if wavelen in normdata[:,0]:
                self.norm=normdata[normdata[:,0]==wavelen,1:3][0]
            else:
                self.norm=self.Normalization()
                norm=array([wavelen,self.norm[0],self.norm[1]])
                norm=norm.reshape((1,3))
                normdata=np.append(normdata, norm, 0)
                np.save(norm_file, normdata)
        except IOError:
            self.norm=self.Normalization()
            normdata=array([wavelen,self.norm[0],self.norm[1]])
            normdata=normdata.reshape((1,3))
            np.save(norm_file, normdata)

        # Calculate the zeroth, first, and second moments of eta' used in the
        # cumulant approximation of the true PSF. Do this separately for the
        # sub- and supercritical regions. This calculation is only needed when
        # the class is initiated.

        def Integrand0_sub(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))
            return integrand

        self.b0sub_norm=qromb(Integrand0_sub,0.0,self.etapmed,1e-4)
        self.b0sub_mean=qromb(lambda etap: Integrand0_sub(etap)*etap,0.0,self.etapmed,1e-4)\
                         /self.b0sub_norm
        self.b0sub_var=qromb(lambda etap: Integrand0_sub(etap)*etap**2,0.0,self.etapmed,1e-4)\
                        /self.b0sub_norm-self.b0sub_mean**2

        def Integrand0_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))
            return integrand

        self.b0real_norm=qromb(Integrand0_real,self.etapmed,self.etapmax,1e-4)
        self.b0real_mean=qromb(lambda etap: Integrand0_real(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                          /self.b0real_norm
        self.b0real_var=qromb(lambda etap: Integrand0_real(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                         /self.b0real_norm-self.b0real_mean**2

        def Integrand0_imag(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))
            return integrand

        self.b0imag_norm=qromb(Integrand0_imag,self.etapmed,self.etapmax,1e-4)
        self.b0imag_mean=qromb(lambda etap: Integrand0_imag(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                          /self.b0imag_norm
        self.b0imag_var=qromb(lambda etap: Integrand0_imag(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                         /self.b0imag_norm-self.b0imag_mean**2

        def Integrand1_sub(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.Epperp(etap)
            return integrand

        self.b1sub_norm=qromb(Integrand1_sub,0.0,self.etapmed,1e-4)
        self.b1sub_mean=qromb(lambda etap: Integrand1_sub(etap)*etap,0.0,self.etapmed,1e-4)\
                         /self.b1sub_norm
        self.b1sub_var=qromb(lambda etap: Integrand1_sub(etap)*etap**2,0.0,self.etapmed,1e-4)\
                        /self.b1sub_norm-self.b1sub_mean**2

        def Integrand1_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.sc5(etap)
            return integrand

        self.b1real_norm=qromb(Integrand1_real,self.etapmed,self.etapmax,1e-4)
        self.b1real_mean=qromb(lambda etap: Integrand1_real(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                          /self.b1real_norm
        self.b1real_var=qromb(lambda etap: Integrand1_real(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                         /self.b1real_norm-self.b1real_mean**2

        def Integrand1_imag(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.sc6(etap)
            return integrand

        self.b1imag_norm=qromb(Integrand1_imag,self.etapmed,self.etapmax,1e-4)
        self.b1imag_mean=qromb(lambda etap: Integrand1_imag(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                          /self.b1imag_norm
        self.b1imag_var=qromb(lambda etap: Integrand1_imag(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                         /self.b1imag_norm-self.b1imag_mean**2

        def Integrand2_sub(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))
            return integrand

        self.b2sub_norm=qromb(Integrand2_sub,0.0,self.etapmed,1e-4)
        self.b2sub_mean=qromb(lambda etap: Integrand2_sub(etap)*etap,0.0,self.etapmed,1e-4)\
                         /self.b2sub_norm
        self.b2sub_var=qromb(lambda etap: Integrand2_sub(etap)*etap**2,0.0,self.etapmed,1e-4)\
                        /self.b2sub_norm-self.b2sub_mean**2

        def Integrand2_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))
            return integrand

        self.b2real_norm=qromb(Integrand2_real,self.etapmed,self.etapmax,1e-4)
        self.b2real_mean=qromb(lambda etap: Integrand2_real(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                          /self.b2real_norm
        self.b2real_var=qromb(lambda etap: Integrand2_real(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                         /self.b2real_norm-self.b2real_mean**2

        def Integrand2_imag(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))
            return integrand

        self.b2imag_norm=qromb(Integrand2_imag,self.etapmed,self.etapmax,1e-4)
        self.b2imag_mean=qromb(lambda etap: Integrand2_imag(etap)*etap,self.etapmed,self.etapmax,1e-4)\
                          /self.b2imag_norm
        self.b2imag_var=qromb(lambda etap: Integrand2_imag(etap)*etap**2,self.etapmed,self.etapmax,1e-4)\
                         /self.b2imag_norm-self.b2imag_mean**2


    # Change and retrieve the radial position value

    def SetRho(self,rho):
        self.rho=rho

    def GetRho(self):
        return self.rho

    # Connect eta and eta0 to etap via Snell's law

    def Eta(self,etap):
        M=self.M
        n_air=self.n_air
        n=self.n
        try:
            eta = np.arcsin(np.clip(M*n_air/n*etap, 0, 1))
        except:
            print(f'{M=}, {n_air=}, {n=}, {etap=}, {M*n_air/n*etap=}')
            raise
        return eta

    def Eta0(self,etap):
        n0=self.n0
        n=self.n
        eta=self.Eta(etap)
        eta0 = np.arccos(sqrt(fabs(1.0-(n/n0*sin(eta))**2)))
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
        n0 = self.n0
        k,k0 = self.k,self.k0
        g=self.gamma(etap)
        value=-(n0*k/k0)*2.0*g**2/(1+g**2)
        return value

    def sc2(self,etap):
        n0=self.n0
        k,k0=self.k,self.k0
        g=self.gamma(etap)
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
        value=(n/n0)*(q/kp)*(k/k0)*2/(1+e**2)
        return value

    def sc6(self,etap):
        n,n0=self.n,self.n0
        k,k0,kp=self.k,self.k0,self.kp
        eta=self.Eta(etap)
        q=k*sin(eta)
        e=self.epsilon(etap)
        value=-(n/n0)*(q/kp)*(k/k0)*2*e/(1+e**2)
        return value

    # The electric field functions

    def Epperp(self,etap):
        n=self.n
        n0=self.n0
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
        kp=self.kp

        eta=self.Eta(etap)
        eta0=self.Eta0(etap)

        eppar=2.0*(n/n0)*(k/kp)*cos(eta)*cos(eta0)*sin(eta)\
               /sin(eta+eta0)/cos(eta0-eta)
        return eppar

    def Espar(self,etap):
        n=self.n
        k=self.k
        k0=self.k0

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

    def Val0_subcrit(self):
        etapmed=self.etapmed


        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))*self.J0(etap)
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)
        return value

    def Val0_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))*self.J0(etap)
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val0_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))*self.J0(etap)
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val1_subcrit(self):
        etapmed=self.etapmed

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.Epperp(etap)*self.J1(etap)
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)
        return value

    def Val1_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.sc5(etap)*self.J1(etap)
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val1_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*self.sc6(etap)*self.J1(etap)
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val2_subcrit(self):
        etapmed=self.etapmed

        def Integrand_real(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))*self.J2(etap)
            return integrand

        value=qromb(Integrand_real,0.0,etapmed,eps=1.0e-4)
        return value

    def Val2_real(self):
        etapmax=self.etapmax
        etapmed=self.etapmed

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))*self.J2(etap)
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    def Val2_imag(self):
        etapmax=self.etapmax
        etapmed=self.etapmed

        def Integrand_sc(etap):
            eta=self.Eta(etap)
            integrand=etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))*self.J2(etap)
            return integrand

        value=qromb(Integrand_sc,etapmed,etapmax,eps=1.0e-4)
        return value

    # Functions for the cumulant expansion

    def Bessel0Cum_subcrit(self):
        rho=self.GetRho()

        norm=self.b0sub_norm
        mom1=self.b0sub_mean
        mom2=self.b0sub_var

        arg=self.kp*rho*mom1

        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish

        j0=jn(0,arg)
        j1=jn(1,arg)
        j2=jn(2,arg)
        value=norm*(j0+0.5*(rho*self.kp)**2*mom2*(j2-j1/arg))
        return value

    def Bessel0Cum_real(self):
        rho=self.GetRho()

        norm=self.b0real_norm
        mom1=self.b0real_mean
        mom2=self.b0real_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        value=norm*(jn(0,arg)+0.5*(rho*self.kp)**2*mom2*(jn(2,arg)-jn(1,arg)/arg))
        return value

    def Bessel0Cum_imag(self):
        rho=self.GetRho()

        norm=self.b0imag_norm
        mom1=self.b0imag_mean
        mom2=self.b0imag_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        value=norm*(jn(0,arg)+0.5*(rho*self.kp)**2*mom2*(jn(2,arg)-jn(1,arg)/arg))
        return value

    def Bessel1Cum_subcrit(self):

        rho=self.GetRho()

        norm=self.b1sub_norm
        mom1=self.b1sub_mean
        mom2=self.b1sub_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel1Cum_real(self):

        rho=self.GetRho()

        norm=self.b1real_norm
        mom1=self.b1real_mean
        mom2=self.b1real_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel1Cum_imag(self):

        rho=self.GetRho()

        norm=self.b1imag_norm
        mom1=self.b1imag_mean
        mom2=self.b1imag_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel2Cum_subcrit(self):

        rho=self.GetRho()

        norm=self.b2sub_norm
        mom1=self.b2sub_mean
        mom2=self.b2sub_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        j2=jn(2,arg)
        value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        return value

    def Bessel2Cum_real(self):

        rho=self.GetRho()

        norm=self.b2real_norm
        mom1=self.b2real_mean
        mom2=self.b2real_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        j2=jn(2,arg)
        value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        return value

    def Bessel2Cum_imag(self):

        rho=self.GetRho()

        norm=self.b2imag_norm
        mom1=self.b2imag_mean
        mom2=self.b2imag_var

        arg=self.kp*rho*mom1
        # Handle the case where arg is zero or very small
        if abs(arg) < 1e-10:
            return norm  # j0(0) = 1, other terms vanish
        j2=jn(2,arg)
        value=norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0,arg)-3*jn(1,arg)/arg+6*j2/arg**2))
        return value

    # Functions for calculation of the normalization constants

    def Intensity_norm_par(self,rho):
        self.SetRho(rho)
        par= ((self.Val0_subcrit()+self.Val0_real())**2+self.Val0_imag()**2\
              +(self.Val2_subcrit()+self.Val2_real())**2+self.Val2_imag()**2)/4.0
        return par

    def Intensity_norm_vert(self,rho):
        self.SetRho(rho)
        vert=((self.Val1_subcrit()+self.Val1_real())**2+self.Val1_imag()**2)
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
        norm_par=qromb(self.NormIntegrand_par,0.0,20000,1e-4)
        norm_vert=qromb(self.NormIntegrand_vert,0.0,20000,1e-4)
        return (norm_par,norm_vert)

    # Calculation of the PSF

    def Intensity_exact(self,rho,phip,alpha,beta):

        self.SetRho(rho)

        par= sin(beta)**2/4.0*\
             ((self.Val0_subcrit()+self.Val0_real())**2+self.Val0_imag()**2\
              +(self.Val2_subcrit()+self.Val2_real())**2+self.Val2_imag()**2\
               -2*cos(2*(phip-alpha))*((self.Val0_subcrit()+self.Val0_real())*\
              (self.Val2_subcrit()+self.Val2_real())\
               +self.Val0_imag()*self.Val2_imag()))

        mix=sin(beta)*cos(beta)*cos(phip-alpha)*\
             ((self.Val1_subcrit()+self.Val1_real())*self.Val2_imag()-
              (self.Val2_subcrit()+self.Val2_real())*self.Val1_imag()-
              (self.Val1_subcrit()+self.Val1_real())*self.Val0_imag()+
              (self.Val0_subcrit()+self.Val0_real())*self.Val1_imag())

        vert=cos(beta)**2*\
              ((self.Val1_subcrit()+self.Val1_real())**2+self.Val1_imag()**2)

        value=par+mix+vert

        # Normalization
        value/=(sin(beta)**2*self.norm[0]+cos(beta)**2*self.norm[1])
        return value

    def Intensity_approx(self,rho,phip,alpha,beta):

        self.SetRho(rho)

        par= sin(beta)**2/4.0*\
             ((self.Bessel0Cum_subcrit()+self.Bessel0Cum_real())**2+self.Bessel0Cum_imag()**2\
              +(self.Bessel2Cum_subcrit()+self.Bessel2Cum_real())**2+self.Bessel2Cum_imag()**2\
               -2*cos(2*(phip-alpha))*((self.Bessel0Cum_subcrit()+self.Bessel0Cum_real())*\
              (self.Bessel2Cum_subcrit()+self.Bessel2Cum_real())\
               +self.Bessel0Cum_imag()*self.Bessel2Cum_imag()))

        mix=sin(beta)*cos(beta)*cos(phip-alpha)*\
             ((self.Bessel1Cum_subcrit()+self.Bessel1Cum_real())*self.Bessel2Cum_imag()-
              (self.Bessel2Cum_subcrit()+self.Bessel2Cum_real())*self.Bessel1Cum_imag()-
              (self.Bessel1Cum_subcrit()+self.Bessel1Cum_real())*self.Bessel0Cum_imag()+
              (self.Bessel0Cum_subcrit()+self.Bessel0Cum_real())*self.Bessel1Cum_imag())

        vert=cos(beta)**2*\
              ((self.Bessel1Cum_subcrit()+self.Bessel1Cum_real())**2+self.Bessel1Cum_imag()**2)

        value=par+mix+vert

        # Normalization
        value/=(sin(beta)**2*self.norm[0]+cos(beta)**2*self.norm[1])
        return value

    def PSF_approx(self,x,y,alpha,beta):
        r=sqrt(x**2+y**2)
        phip = np.arctan2(y, x)%(2*np.pi)
        rho=self.M*r
        value=self.M**2*self.Intensity_approx(rho,phip,alpha,beta)
        return value

    def PSF_exact(self,x,y,alpha,beta):
        r=sqrt(x**2+y**2)
        phip = np.arctan2(y, x)%(2*np.pi)
        rho=self.M*r
        value=self.M**2*self.Intensity_exact(rho,phip,alpha,beta)
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
        #phi,theta,mux_nm,muy_nm,n_photons=x

        if abs(x[0]) > 1:
          print('Problematic x[0]')
          print(x[0])
          print()
        if abs(x[1]) > 1:
          print('Problematic x[1]')
          print(x[1])
          print()
        if abs(x[2]) > 1:
          print('Problematic x[2]')
          print(x[2])
          print()

        phi_conv = np.arctan2(x[1], x[0]) % (2*np.pi)  # Reconvert x[0], x[1] to phi check the order of 0 and 1 needs to be the same as est
        theta_conv = 1.0 * np.arccos(x[2]) % (np.pi/2)  # Reconvert x[2] to theta
        
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
    Estimates the centre coordinates (x and y) and the orientation (theta and phi)
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
    initpix (array) : Array of length 2 of initial values for the centre pixel (ypixel,xpixel)
    deltapix (int)  : The half width of the array to be analyzed

    Functions:
    -------------------------------------------------------------------------
    Estimate (array) : Takes a full pixel array and uses MLEwT to return an
    array of estimates for x,y,b,N,theta,phi where b is the number of background
    photons per pixel and N is the photon number.

    Kim I. Mortensen
    """

    def __init__(self, initvals, initpix, deltapix, psf_generator):

        # Store user input
        
        self.initvals=initvals # Initial values to try in optimisation
        self.initpix=initpix # Central pixe laround which to consider a patch
        self.deltapix = deltapix # Half patch width to operate on

        self.psf_generator = psf_generator
        self._observer = []
        self._verbose = True

    def set_verbose(self, value):
        self._verbose = bool(value)

    def set_observer(self, value):
        assert hasattr(value, 'append')
        self._observer = value

    def Estimate(self, datamatrix):
        
        # Validate initpix
        if not isinstance(self.initpix, (tuple, list)) or len(self.initpix) != 2:
            raise ValueError(f"Invalid initpix: {self.initpix}. Must be a tuple or list of length 2.")
    

        ypix=self.initpix[0]
        xpix=self.initpix[1]
        deltapix=self.deltapix

        #print(f"The values of x, y, delta are 1: ypix={ypix}, xpix={xpix}, deltapix={deltapix}")

        # Ensure ypix and xpix are integers
        if not isinstance(ypix, int) or not isinstance(xpix, int):
            raise TypeError(f"ypix and xpix must be integers. Got ypix={ypix}, xpix={xpix}.")

        # Entire data matrix
        counts=datamatrix
        #print(f"counts inside Estimate() 1: {counts}, shape: {counts.shape}")
 
        # Transformation of initial values
        if not isinstance(self.initvals, (list, np.ndarray)):
            raise ValueError(f"Expected initvals to be a list or np.ndarray of 6 elements, but got {self.initvals} {type(self.initvals)}")
        
        # Reparametrise phi -> sin(phi), cos(phi) this adds a parameter to the fitting
                       #theta -> cos(2*theta)
#        pinit=zeros(6)
#        pinit[0] = np.sin(self.initvals[1]) * np.cos(self.initvals[0])
#        pinit[1] = np.sin(self.initvals[1]) * np.sin(self.initvals[0])
#        pinit[2] = np.cos(self.initvals[1])   #reparametrised theta in cos(2theta)
#        pinit[3] = self.initvals[2]                      # mux_nm
#        pinit[4] = self.initvals[3]                      # muy_nm
#        pinit[5] = self.initvals[4]                      #sqrt(self.initvals[4]) # n_photons
        pinit=zeros(6)
        pinit[0] = self.initvals[0]                      # mux_nm
        pinit[1] = self.initvals[1]                      # mux_nm
        pinit[2] = self.initvals[2]                      # mux_nm
        pinit[3] = self.initvals[3]                      # mux_nm
        pinit[4] = self.initvals[4]                      # mux_nm
        pinit[5] = self.initvals[5]                      # mux_nm

#        print(f"    iinitial theta = {round(np.arccos(pinit[2])*180/np.pi % 90)}")
#        print(f"    initial phi = {round(np.arctan2(pinit[1], pinit[0])*180/np.pi % 360)}")
#        print(f"    initial params = ({round(pinit[0], 4)}, {round(pinit[1], 4)}, {round(pinit[2], 4)})")        



        # Create instance of LogLikelihood object
        ll=LogLikelihood(counts, self.psf_generator)

        # Perform maximization of the log-likelihood using Powell's method
        #bounds = [(0, 2*np.pi), (0, np.pi/2), (-433, 433), (-433, 433), (1e-6, None)]
        bounds = [(-1, 1), (-1, 1), (0, 1), (self.initvals[3] - 200, self.initvals[3] + 200), (self.initvals[4] - 200, self.initvals[4] + 200), (1, 1e10)]

        def constraint(x):
            return x[0]**2 + x[1]**2 + x[2]**2 - 1

        constraints = (
            {'type': 'eq', 'fun': constraint},  # 'eq' means equality constraint
        )

        def zero_hessian(x):
            n = len(x)
            return np.zeros((n, n))

        options = {
            'xtol': 1e-3,
            'ftol': 1e-3,
#            'barrier_tol': 1e-3,
#            'gtol': 1e-1,         # Looser gradient tolerance (default is often 1e-8)
#            'xtol': 1e-1,         # Looser step size tolerance
#            'barrier_tol': 1e-1,  # Looser barrier parameter tolerance
#            'initial_constr_penalty': 1.0,  # Initial constraint penalty
#            'maxiter': 500       # Maximum iterations if needed
        }


        # Store all parameters and scores
        params = []
        scores = []
        max_attempts = 5
        
        # First attempt
        print('    Attempt 1')
        res = minimize(ll.Value, pinit, method='Powell', bounds=bounds, options=options)
#        res = minimize(ll.Value, pinit, method='trust-constr', bounds=bounds, options=options, constraints=constraints)
        first_attempt_params = res.x
        first_attempt_score = ll.Value(first_attempt_params)
        params.append(first_attempt_params)
        scores.append(first_attempt_score)
        
        # Check if theta estimate is near 90 degrees
        phi_est = np.arctan2(first_attempt_params[1], first_attempt_params[0]) % (2*np.pi)
        theta_est = np.arccos(first_attempt_params[2]) % (np.pi/2)
        
        # Hinterer-style repetitions
        attempt = 1
        while attempt < max_attempts:
            conditionA = abs(theta_est - np.pi/2) >= 1*np.pi/180
#            conditionB = abs(phi_est - np.pi) >= 10*np.pi/180
            conditionB = min(abs(phi_est - angle) for angle in [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]) >= 20*np.pi/180
            conditionC = abs(theta_est) < 0
            if conditionA and conditionB and conditionC:
                # Condition satisfied, exit loop
                break
            
            # Increment attempt counter
            attempt += 1
            print(f'    Attempt {attempt}')
            
            # For the second attempt, use the specific modification
            if attempt == 9992:#2:
                next_attempt_pinit = first_attempt_params.copy()
                next_attempt_pinit[0] = np.sin(theta_est)*np.cos(phi_est - np.pi)
                next_attempt_pinit[1] = np.sin(theta_est)*np.sin(phi_est - np.pi)
                next_attempt_pinit[2] = np.cos(theta_est)
            elif attempt == 9993:#3:
                next_attempt_pinit = first_attempt_params.copy()
                next_attempt_pinit[0] = np.sin(theta_est - np.pi/2)*np.cos(phi_est)
                next_attempt_pinit[1] = np.sin(theta_est - np.pi/2)*np.sin(phi_est)
                next_attempt_pinit[2] = np.cos(theta_est - np.pi/2)
            elif attempt == 9994:#4:
                next_attempt_pinit = first_attempt_params.copy()
                next_attempt_pinit[0] = np.sin(theta_est - np.pi/2)*np.cos(phi_est - np.pi)
                next_attempt_pinit[1] = np.sin(theta_est - np.pi/2)*np.sin(phi_est - np.pi)
                next_attempt_pinit[2] = np.cos(theta_est - np.pi/2)
            else:
                # For subsequent attempts, use random parameters
                rand_theta = np.random.uniform(0, np.pi/2)
                rand_phi = np.random.uniform(0, 2 * np.pi)
                next_attempt_pinit = first_attempt_params.copy()
                next_attempt_pinit[0] = np.sin(rand_theta)*np.cos(rand_phi)
                next_attempt_pinit[1] = np.sin(rand_theta)*np.sin(rand_phi)
                next_attempt_pinit[2] = np.cos(rand_theta)
                next_attempt_pinit[3] = np.random.uniform(-51/2, 51/2)
                next_attempt_pinit[4] = np.random.uniform(-51/2, 51/2)
                next_attempt_pinit[5] = np.random.uniform(1, 1e10)
            
            # Perform optimization
            res = minimize(ll.Value, next_attempt_pinit, method='Powell', bounds=bounds, options=options)
#            res = minimize(ll.Value, pinit, method='trust-constr', bounds=bounds, options=options, constraints=constraints)
            next_attempt_params = res.x
            next_attempt_score = ll.Value(next_attempt_params)
            
            # Store results
            params.append(next_attempt_params)
            scores.append(next_attempt_score)
            
            # Update phi_est for next iteration check
            theta_est = np.arccos(next_attempt_params[2]) % (np.pi/2)
            phi_est = np.arctan2(next_attempt_params[1], next_attempt_params[0]) % (2*np.pi)

#        # Print final status
#        if attempt == max_attempts:
#            print(f'    Maximum attempts ({max_attempts}) reached')
#        else:
#            print(f'    Condition satisfied after {attempt} attempts')
        
        # Find the best solution from all attempts
        min_score_index = scores.index(min(scores))
        winning_params = params[min_score_index]
        
        
        
        
        
#        Loads of random additional attempts
#
#        conditionA = abs(theta_est - np.pi/2) <= 0.1*np.pi/180
#        conditionB = any(abs(phi_est - angle) <= 20*np.pi/180 for angle in np.array([0, 90, 180, 270, 360]) * np.pi/180)
#
#        # If it is, repeat until it's not (or max reached)
#        attempt_count = 1
#        max_attempts = 10
#        
#        while (conditionA or conditionB) and attempt_count < max_attempts:
#            attempt_count += 1
#            print(f'    Attempt {attempt_count}')
#            
#            # Generate new initial values
#            random_init = pinit.copy()
#            
#            # First, try using the previous phi but with a different theta
#            # That often works
#            if attempt_count == 2:
#                random_init[0] = np.sin(theta_est)*np.cos(phi_est - np.pi)
#                random_init[1] = np.sin(theta_est)*np.sin(phi_est - np.pi)
#                random_init[2] = np.cos(theta_est)
#            else:
#                # For additional attempts, use random initial values
#                for i in range(len(pinit)):
#                  lower, upper = bounds[i]
#                  random_init[i] = np.random.uniform(lower, upper)
#            
#            res = minimize(ll.Value, random_init, method='Powell', bounds=bounds, options=options)
#            current_params = res.x
#            current_score = ll.Value(current_params)
#            params.append(current_params)
#            scores.append(current_score)
#            
#            # Check if theta estimate is close to 90 degrees
#            phi_est = np.arctan2(current_params[1], current_params[0]) % (2*np.pi)
#            theta_est = 1.0 * np.arccos(current_params[2]) % (np.pi/2)
#            conditionA = abs(theta_est - np.pi/2) <= 0.1*np.pi/180
#            conditionB = any(abs(phi_est - angle) <= 20*np.pi/180 for angle in np.array([0, 90, 180, 270, 360]) * np.pi/180)
#            
#            # If not near 90, finish
#            if not (conditionA or conditionB):
#                break
#        
#        # Find the best solution from all attempts
#        min_score_index = scores.index(min(scores))
#        winning_params = params[min_score_index]



        # Store position estimates relative to initial pixel
        # retransform the parametrised phi and theta in the original ones in the est
        self.phi_conv = np.arctan2(winning_params[1], winning_params[0]) % (2*np.pi)
        self.theta_conv = np.arccos(winning_params[2]) % (np.pi/2)

#        print('theta =', round(self.theta_conv*180/np.pi))
#        print('phi =', round(self.phi_conv*180/np.pi))

        #est[0] = phi_est      # phi reconverted
        #est[1] = theta_est    # theta reconverted
        self.mux_nm = winning_params[3]    # mux_nm
        self.muy_nm = winning_params[4]    # muy_nm
        n_photons = winning_params[5]    # n_photons

        results = [self.phi_conv, self.theta_conv, self.mux_nm, self.muy_nm, n_photons]
        return results
        
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
    per pixel, the centre coordinates (array([x,y])), the azimuthal angle (theta), and the
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
        
        # Matrix of distances from PSF centre to pixel centres
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


if __name__=='__main__':
    home = os.getenv('HOME')
    norm_file = os.path.join(home, 'dipolenorm.npy')
    #example(norm_file)

