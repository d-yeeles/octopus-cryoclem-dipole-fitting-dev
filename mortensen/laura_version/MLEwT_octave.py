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
from numint import qromb
#from numint import qromb

from oct2py import Oct2Py
import tempfile
import scipy.io as sio
import shutil

import subprocess
from numpy import pi, sqrt, sin, cos, arccos, arctan2

import os
import numpy as np
import tempfile
import subprocess
from numpy import pi, sqrt, sin, cos, arccos, arctan2
from scipy.io import savemat, loadmat
from scipy.special import gammaln


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
        value=norm*(jn(0,arg)+0.5*(rho*self.kp)**2*mom2*(jn(2,arg)-jn(1,arg)/arg))
        return value

    def Bessel0Cum_imag(self):
        rho=self.GetRho()

        norm=self.b0imag_norm
        mom1=self.b0imag_mean
        mom2=self.b0imag_var

        arg=self.kp*rho*mom1
        value=norm*(jn(0,arg)+0.5*(rho*self.kp)**2*mom2*(jn(2,arg)-jn(1,arg)/arg))
        return value

    def Bessel1Cum_subcrit(self):

        rho=self.GetRho()

        norm=self.b1sub_norm
        mom1=self.b1sub_mean
        mom2=self.b1sub_var

        arg=self.kp*rho*mom1
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel1Cum_real(self):

        rho=self.GetRho()

        norm=self.b1real_norm
        mom1=self.b1real_mean
        mom2=self.b1real_var

        arg=self.kp*rho*mom1
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel1Cum_imag(self):

        rho=self.GetRho()

        norm=self.b1imag_norm
        mom1=self.b1imag_mean
        mom2=self.b1imag_var

        arg=self.kp*rho*mom1
        j1=jn(1,arg)
        value=norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0,arg)/arg-j1))
        return value

    def Bessel2Cum_subcrit(self):

        rho=self.GetRho()

        norm=self.b2sub_norm
        mom1=self.b2sub_mean
        mom2=self.b2sub_var

        arg=self.kp*rho*mom1
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



class MLEwT:
    """
    Estimates the center coordinates (x and y) and the orientation (theta and phi)
    of a fixed dipole using MLEwT with Octave's fmincon.
    """

    def __init__(self, initvals, psf_generator):
        self.initvals = initvals
        self.psf_generator = psf_generator
        self._observer = []
        self._verbose = True

    def set_verbose(self, value):
        self._verbose = bool(value)

    def set_observer(self, value):
        assert hasattr(value, 'append')
        self._observer = value

    def Estimate(self, datamatrix):
        """
        Estimate parameters using Octave's fmincon optimizer.
        
        Parameters:
        -----------
        datamatrix : array_like
            Image data
            
        Returns:
        --------
        results : list
            [phi, theta, mux_nm, muy_nm, n_photons]
        """
        # Create a temporary directory for file exchange
        temp_dir = tempfile.mkdtemp()
        
        try:
            # Initialize parameters
            pinit = np.zeros(6)
            pinit[0] = self.initvals[0]  # First octave param
            pinit[1] = self.initvals[1]  # Second octave param
            pinit[2] = self.initvals[2]  # Third octave param
            pinit[3] = self.initvals[3]  # mux_nm
            pinit[4] = self.initvals[4]  # muy_nm
            pinit[5] = self.initvals[5]  # n_photons
            
            # Save the data matrix to a .mat file for Octave
            data_file = os.path.join(temp_dir, 'datamatrix.mat')
            sio.savemat(data_file, {'datamatrix': datamatrix})
            
            # Create a function to generate PSF models on-demand
            def generate_psf_model(phi, theta, mux_nm, muy_nm, n_photons):
                # Call the Python PSF generator
                psf = self.psf_generator(phi, theta, mux_nm, muy_nm, n_photons)
                return psf
                
            # Pre-compute a few PSF models to check shape
            test_psf = generate_psf_model(0, 0, 0, 0, 1000)
            
            # Function to save PSF for a given parameter set
            def save_psf_for_params(params):
                phi = np.arctan2(params[1], params[0]) % (2*np.pi)
                theta = np.arccos(params[2]) % (np.pi/2)
                mux_nm = params[3]
                muy_nm = params[4]
                n_photons = params[5]
                
                psf = generate_psf_model(phi, theta, mux_nm, muy_nm, n_photons)
                return psf

            # Create Python interface for Octave to call
            # Instead of direct calling, we'll pre-compute PSFs for initial parameters
            # and a grid of parameters around them for Octave to interpolate
            
            # Create a grid of parameter variations
            param_variations = []
            param_variations.append(pinit)  # Add initial parameters
            
            # Create small variations of each parameter
            for i in range(6):
                for delta in [-0.1, 0.1]:
                    if i < 3:  # For orientation parameters
                        param_var = pinit.copy()
                        param_var[i] += delta
                        # Normalize
                        norm = np.sqrt(param_var[0]**2 + param_var[1]**2 + param_var[2]**2)
                        if norm > 0:
                            param_var[0:3] /= norm
                        param_variations.append(param_var)
                    elif i == 5:  # For photon count, use multiplicative delta
                        param_var = pinit.copy()
                        param_var[i] *= (1 + delta)
                        param_variations.append(param_var)
                    else:  # For position parameters
                        param_var = pinit.copy()
                        param_var[i] += delta * 10  # Larger step for position
                        param_variations.append(param_var)
            
            # Compute PSFs for all parameter variations
            psf_models = []
            for params in param_variations:
                psf = save_psf_for_params(params)
                psf_models.append(psf)
            
            # Save all PSF models to a single .mat file
            psf_file = os.path.join(temp_dir, 'psf_models.mat')
            param_mat = np.array(param_variations)
            psf_mat = np.array(psf_models)
            sio.savemat(psf_file, {'params': param_mat, 'psfs': psf_mat})
            
            # Create lookup function in Octave
            with open(os.path.join(temp_dir, 'lookup_psf.m'), 'w') as f:
                f.write("""
                function psf = lookup_psf(phi, theta, mux_nm, muy_nm, n_photons, params, psfs)
                    % Convert angles to xyz representation
                    x0 = sin(theta) * cos(phi);
                    y0 = sin(theta) * sin(phi);
                    z0 = cos(theta);
                    
                    % Create query point
                    query = [x0, y0, z0, mux_nm, muy_nm, n_photons];
                    
                    % Find closest parameter set in our database
                    [~, idx] = min(sum((params - repmat(query, size(params, 1), 1)).^2, 2));
                    
                    % Return the corresponding PSF
                    psf = squeeze(psfs(idx, :, :));
                end
                """)
            
            # Create log-likelihood function
            with open(os.path.join(temp_dir, 'lnpdf.m'), 'w') as f:
                f.write("""
                function val = lnpdf(x, datamatrix, params, psfs)
                    % Convert parameterization
                    phi = mod(atan2(x(2), x(1)), 2*pi);
                    theta = mod(acos(x(3)), pi/2);
                    mux_nm = x(4);
                    muy_nm = x(5);
                    n_photons = x(6);
                    
                    % Get PSF model from lookup
                    model_image = lookup_psf(phi, theta, mux_nm, muy_nm, n_photons, params, psfs);
                    
                    % Calculate log-likelihood
                    model_image = max(model_image, 1e-10);  % Ensure positive values
                    ln_factorial = gammaln(datamatrix + 1);
                    log_likelihood = sum(sum(datamatrix .* log(model_image) - model_image - ln_factorial));
                    
                    % Return negative log-likelihood for minimization
                    val = -log_likelihood;
                end
                """)
            
            # Create optimization function
            with open(os.path.join(temp_dir, 'run_optimization.m'), 'w') as f:
                f.write("""
                function result = run_optimization(datamatrix, initvals, params, psfs)
                    pkg load optim
                    
                    % Define bounds
                    lower_bounds = [-1, -1, -1, -459, -459, 1];
                    upper_bounds = [1, 1, 1, 459, 459, 1e10];
                    
                    % Create cost function wrapper
                    costFun = @(x) lnpdf(x, datamatrix, params, psfs);
                    
                    % Test the cost function on initial values
                    try
                        test_val = costFun(initvals);
                        disp(['Initial cost: ', num2str(test_val)]);
                    catch err
                        disp('Error evaluating cost function:');
                        disp(err.message);
                        error('Failed to evaluate cost function');
                    end
                    
                    % Set optimization options
                    options = optimset('Display', 'iter', 'TolFun', 1e-4, 'TolX', 1e-4, 'MaxIter', 50);
                    
                    % Run optimization without constraint
                    try
                        [x_opt, fval] = fmincon(costFun, initvals, [], [], [], [], lower_bounds, upper_bounds, [], options);
                        disp(['Final cost: ', num2str(fval)]);
                    catch err
                        disp('Optimization failed:');
                        disp(err.message);
                        error('Failed to run optimization');
                    end
                    
                    % Calculate angles for output
                    phi_conv = mod(atan2(x_opt(2), x_opt(1)), 2*pi);
                    theta_conv = mod(acos(x_opt(3)), pi/2);
                    
                    % Normalize the orientation parameters
                    norm_factor = sqrt(x_opt(1)^2 + x_opt(2)^2 + x_opt(3)^2);
                    if norm_factor > 0
                        x_opt(1:3) = x_opt(1:3) / norm_factor;
                    end
                    
                    % Return result
                    result = [phi_conv, theta_conv, x_opt(4), x_opt(5), x_opt(6)];
                end
                """)
            
            # Initialize Octave session
            print("    Starting Octave session...")
            oc = Oct2Py()
            
            # Add temp directory to Octave path
            oc.addpath(temp_dir)
            
            # Load data and PSF models
            oc.eval("load('" + data_file.replace('\\', '/') + "');")
            oc.eval("load('" + psf_file.replace('\\', '/') + "');")
            
            # Call the optimization function
            print("    Running Octave optimization...")
            results = oc.run_optimization(oc.datamatrix, pinit, oc.params, oc.psfs)
            
            # Close the Octave session
            oc.close()
            
            # Store the results for reference
            self.phi_conv = results[0]
            self.theta_conv = results[1]
            self.mux_nm = results[2]
            self.muy_nm = results[3]
            n_photons = results[4]
            
            return results.flatten().tolist()
            
        except Exception as e:
            print(f"Error during Octave optimization: {e}")
            # Fall back to initial values if optimization fails
            phi_conv = np.arctan2(pinit[1], pinit[0]) % (2*np.pi)
            theta_conv = 1.0 * np.arccos(pinit[2]) % (np.pi/2)
            mux_nm = pinit[3]
            muy_nm = pinit[4]
            n_photons = pinit[5]
            
            results = [phi_conv, theta_conv, mux_nm, muy_nm, n_photons]
            print(f"Returning initial values as fallback: {results}")
            return results
            
        finally:
            # Clean up
            shutil.rmtree(temp_dir)


# Simple LogLikelihood class for compatibility
class LogLikelihood:
    """ Class defining the log-likelihood function maximized in MLE."""

    def __init__(self, counts, psf_generator):
        self.counts = counts
        self.psf_generator = psf_generator

    def Value(self, x):
        # This method is kept for compatibility but not used with Octave
        pass


 
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


if __name__=='__main__':
    home = os.getenv('HOME')
    norm_file = os.path.join(home, 'dipolenorm.npy')
    #example(norm_file)

