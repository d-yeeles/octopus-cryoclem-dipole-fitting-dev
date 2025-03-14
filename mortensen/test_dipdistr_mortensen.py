"""
Python 3 Compatible Mortensen PSF Model (dipdistr)

This is the original Mortensen PSF model implementation updated to Python 3.

Original: September 8, 2015, Kim I. Mortensen
Updated to Python 3: February 25, 2025
"""

from scipy.special import jn
from numpy import *
import numpy as np
from numeric_utilities import qromb, find_nearest

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

    def __init__(self, wl, n, n0, M, NA):
        self.wl = wl          # Wavelength in nm (emission peak in buffer)
        self.NA = NA          # Numerical aperture
        
        self.M = M            # Magnification
        self.n = n            # Refractive index of glass/oil
        self.np = 1.0         # Refractive index of air in lab (vacuum)
        self.n0 = n0          # Refractive index of sample medium (water)

        self.kp = 2*pi/wl     # Wavevector amplitude in air (vacuum)
        self.k0 = n0*self.kp  # Wavevector amplitude in buffer
        self.k = n*self.kp    # Wavevector amplitude in glass

        self.etapmed = n0/M   # Integration limits
        self.etapmax = NA/M

        # Calculate or load normalization constants
        try:
            self.norm = numpy.loadtxt('dipoletablenorm.dat')
        except:
            self.deltaz = -1.0
            self.norm = self.Normalization()
            numpy.savetxt('dipoletablenorm.dat', self.norm)

        # Tabulate the PSF approximations at this range of focal values 
        self.focusvals = arange(-150, 150, 1.0)
        self.nfoci = len(self.focusvals)
        nf = self.nfoci

        self.b0subreal_norm = zeros(nf)
        self.b0subreal_mean = zeros(nf)
        self.b0subreal_var = zeros(nf)
        self.b0subimag_norm = zeros(nf)
        self.b0subimag_mean = zeros(nf)
        self.b0subimag_var = zeros(nf)
        self.b0real_norm = zeros(nf)
        self.b0real_mean = zeros(nf)
        self.b0real_var = zeros(nf)
        self.b0imag_norm = zeros(nf)
        self.b0imag_mean = zeros(nf)
        self.b0imag_var = zeros(nf)

        self.b1subreal_norm = zeros(nf)
        self.b1subreal_mean = zeros(nf)
        self.b1subreal_var = zeros(nf)
        self.b1subimag_norm = zeros(nf)
        self.b1subimag_mean = zeros(nf)
        self.b1subimag_var = zeros(nf)
        self.b1real_norm = zeros(nf)
        self.b1real_mean = zeros(nf)
        self.b1real_var = zeros(nf)
        self.b1imag_norm = zeros(nf)
        self.b1imag_mean = zeros(nf)
        self.b1imag_var = zeros(nf)

        self.b2subreal_norm = zeros(nf)
        self.b2subreal_mean = zeros(nf)
        self.b2subreal_var = zeros(nf)
        self.b2subimag_norm = zeros(nf)
        self.b2subimag_mean = zeros(nf)
        self.b2subimag_var = zeros(nf)
        self.b2real_norm = zeros(nf)
        self.b2real_mean = zeros(nf)
        self.b2real_var = zeros(nf)
        self.b2imag_norm = zeros(nf)
        self.b2imag_mean = zeros(nf)
        self.b2imag_var = zeros(nf)
        
        
        for nf in range(self.nfoci):
            self.deltaz = self.focusvals[nf]       

            # Calculate the zeroth, first, and second moments of eta' used in the
            # cumulant approximation of the true PSF. Do this separately for the
            # sub- and supercritical regions. This calculation is only needed when
            # the class is initiated.
            
            def Integrand0_subreal(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))\
                           *cos(self.k*self.deltaz*cos(eta))
                return integrand

            self.b0subreal_norm[nf] = qromb(Integrand0_subreal, 0.0, self.etapmed, 1e-4)
            self.b0subreal_mean[nf] = qromb(lambda etap: Integrand0_subreal(etap)*etap, 0.0, self.etapmed, 1e-4)\
                             /self.b0subreal_norm[nf]
            self.b0subreal_var[nf] = qromb(lambda etap: Integrand0_subreal(etap)*etap**2, 0.0, self.etapmed, 1e-4)\
                            /self.b0subreal_norm[nf]-self.b0subreal_mean[nf]**2
           
            def Integrand0_subimag(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            if abs(self.deltaz) >= 1.0:
                self.b0subimag_norm[nf] = qromb(Integrand0_subimag, 0.0, self.etapmed, 1e-4)
                self.b0subimag_mean[nf] = qromb(lambda etap: Integrand0_subimag(etap)*etap, 0.0, self.etapmed, 1e-4)\
                                 /self.b0subimag_norm[nf]
                self.b0subimag_var[nf] = qromb(lambda etap: Integrand0_subimag(etap)*etap**2, 0.0, self.etapmed, 1e-4)\
                                /self.b0subimag_norm[nf]-self.b0subimag_mean[nf]**2

            def Integrand0_real(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           -etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b0real_norm[nf] = qromb(Integrand0_real, self.etapmed, self.etapmax, 1e-4)
            self.b0real_mean[nf] = qromb(lambda etap: Integrand0_real(etap)*etap, self.etapmed, self.etapmax, 1e-4)\
                              /self.b0real_norm[nf]
            self.b0real_var[nf] = qromb(lambda etap: Integrand0_real(etap)*etap**2, self.etapmed, self.etapmax, 1e-4)\
                             /self.b0real_norm[nf]-self.b0real_mean[nf]**2

            def Integrand0_imag(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           +etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b0imag_norm[nf] = qromb(Integrand0_imag, self.etapmed, self.etapmax, 1e-4)
            self.b0imag_mean[nf] = qromb(lambda etap: Integrand0_imag(etap)*etap, self.etapmed, self.etapmax, 1e-4)\
                              /self.b0imag_norm[nf]
            self.b0imag_var[nf] = qromb(lambda etap: Integrand0_imag(etap)*etap**2, self.etapmed, self.etapmax, 1e-4)\
                             /self.b0imag_norm[nf]-self.b0imag_mean[nf]**2

            def Integrand1_subreal(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*self.Epperp(etap)\
                           *cos(self.k*self.deltaz*cos(eta))
                return integrand

            self.b1subreal_norm[nf] = qromb(Integrand1_subreal, 0.0, self.etapmed, 1e-4)
            self.b1subreal_mean[nf] = qromb(lambda etap: Integrand1_subreal(etap)*etap, 0.0, self.etapmed, 1e-4)\
                             /self.b1subreal_norm[nf]
            self.b1subreal_var[nf] = qromb(lambda etap: Integrand1_subreal(etap)*etap**2, 0.0, self.etapmed, 1e-4)\
                            /self.b1subreal_norm[nf]-self.b1subreal_mean[nf]**2

            def Integrand1_subimag(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*self.Epperp(etap)\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            if abs(self.deltaz) >= 1.0:
                self.b1subimag_norm[nf] = qromb(Integrand1_subimag, 0.0, self.etapmed, 1e-4)
                self.b1subimag_mean[nf] = qromb(lambda etap: Integrand1_subimag(etap)*etap, 0.0, self.etapmed, 1e-4)\
                                 /self.b1subimag_norm[nf]
                self.b1subimag_var[nf] = qromb(lambda etap: Integrand1_subimag(etap)*etap**2, 0.0, self.etapmed, 1e-4)\
                                /self.b1subimag_norm[nf]-self.b1subimag_mean[nf]**2

            def Integrand1_real(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*self.sc5(etap)\
                           *cos(self.k*self.deltaz*cos(eta))\
                           -etap/sqrt(cos(eta))*self.sc6(etap)\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b1real_norm[nf] = qromb(Integrand1_real, self.etapmed, self.etapmax, 1e-4)
            self.b1real_mean[nf] = qromb(lambda etap: Integrand1_real(etap)*etap, self.etapmed, self.etapmax, 1e-4)\
                              /self.b1real_norm[nf]
            self.b1real_var[nf] = qromb(lambda etap: Integrand1_real(etap)*etap**2, self.etapmed, self.etapmax, 1e-4)\
                             /self.b1real_norm[nf]-self.b1real_mean[nf]**2

            def Integrand1_imag(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*self.sc6(etap)\
                           *cos(self.k*self.deltaz*cos(eta))\
                           +etap/sqrt(cos(eta))*self.sc5(etap)\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b1imag_norm[nf] = qromb(Integrand1_imag, self.etapmed, self.etapmax, 1e-4)
            self.b1imag_mean[nf] = qromb(lambda etap: Integrand1_imag(etap)*etap, self.etapmed, self.etapmax, 1e-4)\
                              /self.b1imag_norm[nf]
            self.b1imag_var[nf] = qromb(lambda etap: Integrand1_imag(etap)*etap**2, self.etapmed, self.etapmax, 1e-4)\
                             /self.b1imag_norm[nf]-self.b1imag_mean[nf]**2


            def Integrand2_subreal(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))\
                           *cos(self.k*self.deltaz*cos(eta))
                return integrand

            self.b2subreal_norm[nf] = qromb(Integrand2_subreal, 0.0, self.etapmed, 1e-4)
            self.b2subreal_mean[nf] = qromb(lambda etap: Integrand2_subreal(etap)*etap, 0.0, self.etapmed, 1e-4)\
                             /self.b2subreal_norm[nf]
            self.b2subreal_var[nf] = qromb(lambda etap: Integrand2_subreal(etap)*etap**2, 0.0, self.etapmed, 1e-4)\
                            /self.b2subreal_norm[nf]-self.b2subreal_mean[nf]**2


            def Integrand2_subimag(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            if abs(self.deltaz) >= 1.0:
                self.b2subimag_norm[nf] = qromb(Integrand2_subimag, 0.0, self.etapmed, 1e-4)
                self.b2subimag_mean[nf] = qromb(lambda etap: Integrand2_subimag(etap)*etap, 0.0, self.etapmed, 1e-4)\
                                 /self.b2subimag_norm[nf]
                self.b2subimag_var[nf] = qromb(lambda etap: Integrand2_subimag(etap)*etap**2, 0.0, self.etapmed, 1e-4)\
                                /self.b2subimag_norm[nf]-self.b2subimag_mean[nf]**2


            def Integrand2_real(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           -etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b2real_norm[nf] = qromb(Integrand2_real, self.etapmed, self.etapmax, 1e-4)
            self.b2real_mean[nf] = qromb(lambda etap: Integrand2_real(etap)*etap, self.etapmed, self.etapmax, 1e-4)\
                              /self.b2real_norm[nf]
            self.b2real_var[nf] = qromb(lambda etap: Integrand2_real(etap)*etap**2, self.etapmed, self.etapmax, 1e-4)\
                             /self.b2real_norm[nf]-self.b2real_mean[nf]**2

            def Integrand2_imag(etap):
                eta = self.Eta(etap)
                integrand = etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))\
                           *cos(self.k*self.deltaz*cos(eta))\
                           +etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))\
                           *sin(self.k*self.deltaz*cos(eta))
                return integrand

            self.b2imag_norm[nf] = qromb(Integrand2_imag, self.etapmed, self.etapmax, 1e-4)
            self.b2imag_mean[nf] = qromb(lambda etap: Integrand2_imag(etap)*etap, self.etapmed, self.etapmax, 1e-4)\
                              /self.b2imag_norm[nf]
            self.b2imag_var[nf] = qromb(lambda etap: Integrand2_imag(etap)*etap**2, self.etapmed, self.etapmax, 1e-4)\
                             /self.b2imag_norm[nf]-self.b2imag_mean[nf]**2

    # Change and retrieve the focus
    def SetFocus(self, deltaz):
        self.deltaz = deltaz

    def GetFocus(self):
        return self.deltaz
 
    # Change and retrieve the radial position value
    def SetRho(self, rho):
        self.rho = rho

    def GetRho(self):
        return self.rho

    # Connect eta and eta0 to etap via Snell's law
    def Eta(self, etap):
        M = self.M
        np = self.np
        n = self.n
        eta = arcsin(M*np/n*etap)
        return eta

    def Eta0(self, etap):
        M = self.M
        np = self.np
        n0 = self.n0
        n = self.n
        eta = self.Eta(etap)
        eta0 = arccos(sqrt(fabs(1.0-(n/n0*sin(eta))**2)))
        return eta0
    
    # The Fresnel transmission coefficients
    def Ts(self, etap):
        eta = self.Eta(etap)
        eta0 = self.Eta0(etap)
        ts = 2*cos(eta0)*sin(eta)/sin(eta0+eta)
        return ts

    def Tp(self, etap):
        eta = self.Eta(etap)
        eta0 = self.Eta0(etap)
        ts = self.Ts(etap)
        tp = ts/cos(eta0-eta)
        return tp

    # z-components of wavevectors
    def W(self, etap):
        eta = self.Eta(etap)
        k = self.k
        w = cos(eta)*k
        return w

    def W0(self, etap):
        eta0 = self.Eta0(etap)
        k0 = self.k0
        w0 = cos(eta0)*k0
        return w0

    # Support functions
    def com(self, etap):
        n = self.n
        n0 = self.n0
        eta = self.Eta(etap)
        value = sqrt(fabs(1-(n/n0*sin(eta))**2))
        return value

    def gamma(self, etap):
        n = self.n
        n0 = self.n0
        eta = self.Eta(etap)
        c = self.com(etap)
        value = (n/n0)*cos(eta)/c
        return value

    def delta(self, etap):
        n = self.n
        n0 = self.n0
        eta = self.Eta(etap)
        c = self.com(etap)
        value = (n0/n)*cos(eta)/c
        return value

    def epsilon(self, etap):
        n = self.n
        n0 = self.n0
        eta = self.Eta(etap)
        c = self.com(etap)
        value = (n/n0)*c/cos(eta)
        return value

    # Integrands for super-critical angles
    def sc1(self, etap):
        n, n0 = self.n, self.n0
        k, k0, kp = self.k, self.k0, self.kp
        g = self.gamma(etap)
        c = self.com(etap)
        value = -(n0*k/k0)*2.0*g**2/(1+g**2)
        return value

    def sc2(self, etap):
        n, n0 = self.n, self.n0
        k, k0, kp = self.k, self.k0, self.kp
        g = self.gamma(etap)
        c = self.com(etap)
        value = (n0*k/k0)*2.0*g/(1+g**2)
        return value

    def sc3(self, etap):
        n, n0 = self.n, self.n0
        k, kp = self.k, self.kp
        d = self.delta(etap)
        c = self.com(etap)
        value = 2*(n/n0)*(k/kp)*d*c/(1+d**2)
        return value

    def sc4(self, etap):
        n, n0 = self.n, self.n0
        k, kp = self.k, self.kp
        d = self.delta(etap)
        c = self.com(etap)
        value = 2*(n/n0)*(k/kp)*c*d**2/(1+d**2)
        return value

    def sc5(self, etap):
        n, n0 = self.n, self.n0
        k, k0, kp = self.k, self.k0, self.kp
        eta = self.Eta(etap)
        q = k*sin(eta)
        e = self.epsilon(etap)
        c = self.com(etap)
        value = (n/n0)*(q/kp)*(k/k0)*2/(1+e**2)
        return value

def sc6(self, etap):
        n, n0 = self.n, self.n0
        k, k0, kp = self.k, self.k0, self.kp
        eta = self.Eta(etap)
        q = k*sin(eta)
        e = self.epsilon(etap)
        c = self.com(etap)
        value = -(n/n0)*(q/kp)*(k/k0)*2*e/(1+e**2)
        return value

    # The electric field functions
    def Epperp(self, etap):
        n = self.n
        n0 = self.n0
        np = self.np
        k0 = self.k0
        k = self.k
        kp = self.kp

        eta = self.Eta(etap)
        eta0 = self.Eta0(etap)
        q = k*sin(eta)
                        
        epperp = 2.0*(q/kp)*(n/n0)*(k/k0)*sin(eta)*cos(eta)\
                /sin(eta+eta0)/cos(eta0-eta)
        return epperp

    def Eppar(self, etap):
        n = self.n
        n0 = self.n0
        k = self.k
        k0 = self.k0
        kp = self.kp

        eta = self.Eta(etap)
        eta0 = self.Eta0(etap)

        eppar = 2.0*(n/n0)*(k/kp)*cos(eta)*cos(eta0)*sin(eta)\
               /sin(eta+eta0)/cos(eta0-eta)
        return eppar
        
    def Espar(self, etap):
        n = self.n
        n0 = self.n0
        k = self.k
        k0 = self.k0
        kp = self.kp

        eta = self.Eta(etap)
        eta0 = self.Eta0(etap)
    
        espar = -2.0*n*(k/k0)*cos(eta)*sin(eta)/sin(eta+eta0)
        return espar

    # Bessel functions of appropriate arguments
    def J0(self, etap):
        kp = self.kp
        rho = self.GetRho()
        j0 = jn(0, kp*rho*etap)
        return j0

    def J1(self, etap):
        kp = self.kp
        rho = self.GetRho()
        j1 = jn(1, kp*rho*etap)
        return j1
    
    def J2(self, etap):
        kp = self.kp
        rho = self.GetRho()
        j2 = jn(2, kp*rho*etap)
        return j2
    
    # Calculation of the intensity
    def Val0_subreal(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_real(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))*self.J0(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_real, 0.0, etapmed, eps=1.0e-4)            
        return value

    def Val0_subimag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_real(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)-self.Espar(etap))*self.J0(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_real, 0.0, etapmed, eps=1.0e-4)            
        return value

    def Val0_real(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_sc(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))*\
                       self.J0(etap)*cos(self.kp*self.deltaz*cos(eta))\
                       -etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))*\
                       self.J0(etap)*sin(self.kp*self.deltaz*cos(eta))
            
            return integrand

        value = qromb(Integrand_sc, etapmed, etapmax, eps=1.0e-4)
        return value
        
def Val0_imag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_sc(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.sc4(etap)-self.sc2(etap))*\
                       self.J0(etap)*cos(self.kp*self.deltaz*cos(eta))\
                       +etap/sqrt(cos(eta))*(self.sc3(etap)-self.sc1(etap))*\
                       self.J0(etap)*sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_sc, etapmed, etapmax, eps=1.0e-4)
        return value

    def Val1_subreal(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_real(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*self.Epperp(etap)*self.J1(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_real, 0.0, etapmed, eps=1.0e-4)            
        return value

    def Val1_subimag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_real(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*self.Epperp(etap)*self.J1(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_real, 0.0, etapmed, eps=1.0e-4)            
        return value

    def Val1_real(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_sc(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*self.sc5(etap)*self.J1(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       -etap/sqrt(cos(eta))*self.sc6(etap)*self.J1(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_sc, etapmed, etapmax, eps=1.0e-4)
        return value

    def Val1_imag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_sc(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*self.sc6(etap)*self.J1(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       +etap/sqrt(cos(eta))*self.sc5(etap)*self.J1(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_sc, etapmed, etapmax, eps=1.0e-4)
        return value

    def Val2_subreal(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_real(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))*self.J2(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_real, 0.0, etapmed, eps=1.0e-4)            
        return value

    def Val2_subimag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_real(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.Eppar(etap)+self.Espar(etap))*self.J2(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_real, 0.0, etapmed, eps=1.0e-4)            
        return value
        
def Val2_real(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_sc(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))*self.J2(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       -etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))*self.J2(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_sc, etapmed, etapmax, eps=1.0e-4)
        return value

    def Val2_imag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        def Integrand_sc(etap):
            eta = self.Eta(etap)
            integrand = etap/sqrt(cos(eta))*(self.sc4(etap)+self.sc2(etap))*self.J2(etap)\
                       *cos(self.kp*self.deltaz*cos(eta))\
                       +etap/sqrt(cos(eta))*(self.sc3(etap)+self.sc1(etap))*self.J2(etap)\
                       *sin(self.kp*self.deltaz*cos(eta))
            return integrand

        value = qromb(Integrand_sc, etapmed, etapmax, eps=1.0e-4)
        return value

    # Functions for the cumulant expansion
    def Bessel0Cum_subreal(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b0subreal_norm[self.focusindex]
        mom1 = self.b0subreal_mean[self.focusindex]
        mom2 = self.b0subreal_var[self.focusindex]

        arg = self.kp*rho*mom1
        j0 = jn(0, arg)
        j1 = jn(1, arg)
        j2 = jn(2, arg)
        value = norm*(j0+0.5*(rho*self.kp)**2*mom2*(j2-j1/arg))
        return value

    def Bessel0Cum_subimag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b0subimag_norm[self.focusindex]
        mom1 = self.b0subimag_mean[self.focusindex]
        mom2 = self.b0subimag_var[self.focusindex]

        if mom1 != 0.0:
            arg = self.kp*rho*mom1
            j0 = jn(0, arg)
            j1 = jn(1, arg)
            j2 = jn(2, arg)
            value = norm*(j0+0.5*(rho*self.kp)**2*mom2*(j2-j1/arg))
        else:
            value = 0.0
        return value

    def Bessel0Cum_real(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b0real_norm[self.focusindex]
        mom1 = self.b0real_mean[self.focusindex]
        mom2 = self.b0real_var[self.focusindex]

        arg = self.kp*rho*mom1
        value = norm*(jn(0, arg)+0.5*(rho*self.kp)**2*mom2*(jn(2, arg)-jn(1, arg)/arg))
        return value

    def Bessel0Cum_imag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b0imag_norm[self.focusindex]
        mom1 = self.b0imag_mean[self.focusindex]
        mom2 = self.b0imag_var[self.focusindex]

        arg = self.kp*rho*mom1
        value = norm*(jn(0, arg)+0.5*(rho*self.kp)**2*mom2*(jn(2, arg)-jn(1, arg)/arg))
        return value
        
def Bessel1Cum_subreal(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b1subreal_norm[self.focusindex]
        mom1 = self.b1subreal_mean[self.focusindex]
        mom2 = self.b1subreal_var[self.focusindex]

        arg = self.kp*rho*mom1
        j1 = jn(1, arg)
        value = norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0, arg)/arg-j1))
        return value

    def Bessel1Cum_subimag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b1subimag_norm[self.focusindex]
        mom1 = self.b1subimag_mean[self.focusindex]
        mom2 = self.b1subimag_var[self.focusindex]

        if mom1 != 0.0:
            arg = self.kp*rho*mom1
            j1 = jn(1, arg)
            value = norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0, arg)/arg-j1))
        else:
            value = 0.0
        return value

    def Bessel1Cum_real(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b1real_norm[self.focusindex]
        mom1 = self.b1real_mean[self.focusindex]
        mom2 = self.b1real_var[self.focusindex]

        arg = self.kp*rho*mom1
        j1 = jn(1, arg)
        value = norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0, arg)/arg-j1))
        return value

    def Bessel1Cum_imag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b1imag_norm[self.focusindex]
        mom1 = self.b1imag_mean[self.focusindex]
        mom2 = self.b1imag_var[self.focusindex]

        arg = self.kp*rho*mom1
        j1 = jn(1, arg)
        value = norm*(j1+0.5*(rho*self.kp)**2*mom2*(2*j1/arg**2-jn(0, arg)/arg-j1))
        return value

    def Bessel2Cum_subreal(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b2subreal_norm[self.focusindex]
        mom1 = self.b2subreal_mean[self.focusindex]
        mom2 = self.b2subreal_var[self.focusindex]

        arg = self.kp*rho*mom1
        j2 = jn(2, arg)
        value = norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0, arg)-3*jn(1, arg)/arg+6*j2/arg**2))
        return value

    def Bessel2Cum_subimag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b2subimag_norm[self.focusindex]
        mom1 = self.b2subimag_mean[self.focusindex]
        mom2 = self.b2subimag_var[self.focusindex]

        if mom1 != 0.0:
            arg = self.kp*rho*mom1
            j2 = jn(2, arg)
            value = norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                        (jn(0, arg)-3*jn(1, arg)/arg+6*j2/arg**2))
        else:
            value = 0.0
        return value

    def Bessel2Cum_real(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        kp = self.kp
        rho = self.GetRho()

        norm = self.b2real_norm[self.focusindex]
        mom1 = self.b2real_mean[self.focusindex]
        mom2 = self.b2real_var[self.focusindex]

        arg = self.kp*rho*mom1
        j2 = jn(2, arg)
        value = norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0, arg)-3*jn(1, arg)/arg+6*j2/arg**2))
        return value

    def Bessel2Cum_imag(self):
        etapmax = self.etapmax
        etapmed = self.etapmed
        M = self.M
        np = self.np
        n = self.n

        rho = self.GetRho()

        norm = self.b2imag_norm[self.focusindex]
        mom1 = self.b2imag_mean[self.focusindex]
        mom2 = self.b2imag_var[self.focusindex]

        arg = self.kp*rho*mom1
        j2 = jn(2, arg)
        value = norm*(j2+0.5*(rho*self.kp)**2*mom2*\
                    (jn(0, arg)-3*jn(1, arg)/arg+6*j2/arg**2))
        return value
        
# Functions for calculation of the normalization constants
    def Intensity_norm_par(self, rho):
        self.SetRho(rho)
        par = ((self.Val0_subreal()+self.Val0_real())**2+(self.Val0_subimag()+self.Val0_imag())**2\
              +(self.Val2_subreal()+self.Val2_real())**2+(self.Val2_subimag()+self.Val2_imag())**2)\
              /4.0
        return par

    def Intensity_norm_vert(self, rho):
        self.SetRho(rho)
        vert = ((self.Val1_subreal()+self.Val1_real())**2+(self.Val1_subimag()+self.Val1_imag())**2)
        return vert

    def NormIntegrand_par(self, r):
        rho = self.M*r
        value = zeros(len(rho))
        for i in range(len(rho)):
            value[i] = r[i]*self.Intensity_norm_par(rho[i])
        value *= (2*pi*self.M**2)
        return value

    def NormIntegrand_vert(self, r):
        rho = self.M*r
        value = zeros(len(rho))
        for i in range(len(rho)):
            value[i] = r[i]*self.Intensity_norm_vert(rho[i])
        value *= (2*pi*self.M**2)
        return value

    def Normalization(self):
        print("Calculating normalization constants...")
        norm_par = qromb(self.NormIntegrand_par, 0.0, 20000, 1e-4)
        norm_vert = qromb(self.NormIntegrand_vert, 0.0, 20000, 1e-4)             
        return (norm_par, norm_vert)
    
    # Calculation of the PSF
    def Intensity_exact(self, rho, phip, alpha, beta, deltaz):
        self.focusindex = find_nearest(self.focusvals, deltaz)
        self.deltaz = deltaz
        
        self.SetRho(rho)

        par = sin(beta)**2/4.0*\
             ((self.Val0_subreal()+self.Val0_real())**2+(self.Val0_subimag()+self.Val0_imag())**2\
              +(self.Val2_subreal()+self.Val2_real())**2+(self.Val2_subimag()+self.Val2_imag())**2\
               -2*cos(2*(phip-alpha))*((self.Val0_subreal()+self.Val0_real())*\
              (self.Val2_subreal()+self.Val2_real())\
               +(self.Val0_subimag()+self.Val0_imag())*(self.Val2_subimag()+self.Val2_imag())))

        mix = sin(beta)*cos(beta)*cos(phip-alpha)*\
             ((self.Val1_subreal()+self.Val1_real())*(self.Val2_subimag()+self.Val2_imag())-
              (self.Val2_subreal()+self.Val2_real())*(self.Val1_subimag()+self.Val1_imag())-
              (self.Val1_subreal()+self.Val1_real())*(self.Val0_subimag()+self.Val0_imag())+
              (self.Val0_subreal()+self.Val0_real())*(self.Val1_subimag()+self.Val1_imag()))

        vert = cos(beta)**2*\
              ((self.Val1_subreal()+self.Val1_real())**2+(self.Val1_subimag()+self.Val1_imag())**2)
        
        value = par+mix+vert

        # Normalization
        value /= (sin(beta)**2*self.norm[0]+cos(beta)**2*self.norm[1])
        return value

def Intensity_approx(self, rho, phip, alpha, beta, deltaz):
        self.focusindex = find_nearest(self.focusvals, deltaz)
        self.deltaz = deltaz
        
        self.SetRho(rho)

        par = sin(beta)**2/4.0*\
             ((self.Bessel0Cum_subreal()+self.Bessel0Cum_real())**2+\
              (self.Bessel0Cum_subimag()+ self.Bessel0Cum_imag())**2\
              +(self.Bessel2Cum_subreal()+self.Bessel2Cum_real())**2+\
              (self.Bessel2Cum_subimag()+self.Bessel2Cum_imag())**2\
               -2*cos(2*(phip-alpha))*((self.Bessel0Cum_subreal()+self.Bessel0Cum_real())*\
              (self.Bessel2Cum_subreal()+self.Bessel2Cum_real())\
               +(self.Bessel0Cum_subimag()+self.Bessel0Cum_imag())*\
                                       (self.Bessel2Cum_subimag()+self.Bessel2Cum_imag())))

        mix = sin(beta)*cos(beta)*cos(phip-alpha)*\
             ((self.Bessel1Cum_subreal()+self.Bessel1Cum_real())*\
              (self.Bessel2Cum_subimag()+self.Bessel2Cum_imag())-\
              (self.Bessel2Cum_subreal()+self.Bessel2Cum_real())*\
              (self.Bessel1Cum_subimag()+self.Bessel1Cum_imag())-\
              (self.Bessel1Cum_subreal()+self.Bessel1Cum_real())*\
              (self.Bessel0Cum_subimag()+self.Bessel0Cum_imag())+\
              (self.Bessel0Cum_subreal()+self.Bessel0Cum_real())*\
              (self.Bessel1Cum_subimag()+self.Bessel1Cum_imag()))

        vert = cos(beta)**2*\
              ((self.Bessel1Cum_subreal()+self.Bessel1Cum_real())**2+\
               (self.Bessel1Cum_subimag()+self.Bessel1Cum_imag())**2)
        
        value = par+mix+vert

        # Normalization
        value /= (sin(beta)**2*self.norm[0]+cos(beta)**2*self.norm[1])
        return value

    def PSF_approx(self, x, y, alpha, beta, deltaz):
        # Transform alpha to define optical axis in direction from camera towards objective
        alpha += pi
        self.deltaz = deltaz
        r = sqrt(x**2+y**2)
        if (x < 0.0): 
            phip = pi-arctan(-y/x)
        else: 
            phip = arctan(y/x)
        rho = self.M*r
        value = self.M**2*self.Intensity_approx(rho, phip, alpha, beta, deltaz)
        return value

    def PSF_exact(self, x, y, alpha, beta, deltaz):
        # Transform alpha to define optical axis in direction from camera towards objective
        alpha += pi
        self.deltaz = deltaz
        r = sqrt(x**2+y**2)
        if (x < 0.0): 
            phip = pi-arctan(-y/x)
        else: 
            phip = arctan(y/x)
        rho = self.M*r
        value = self.M**2*self.Intensity_exact(rho, phip, alpha, beta, deltaz)
        return value
        

