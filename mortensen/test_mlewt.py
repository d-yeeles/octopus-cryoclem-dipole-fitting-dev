"""
Maximum Likelihood Estimation with Model Selection Support.

This module provides the MLEwT class with support for both Mortensen and Hinterer PSF models.

Updated to Python 3: February 25, 2025
"""

from scipy.special import i1
from scipy.optimize import fmin_powell
from numpy import *
import numpy as np
from numeric_utilities import gauss, erf

# Import PSF models
from dipdistr_mortensen import dipdistr
from dipdistr_hinterer import dipdistr_hinterer

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
        model (str)     : PSF model to use ('mortensen' or 'hinterer')

        Functions:
        ------------------------------------------------------------
        Value           : Calculates the negative log-likelihood
    """

    def __init__(self, counts, pw, wl, NA, n, n0, M, pinit, Sfloor, alpha, sigma, model='mortensen'):

        self.counts = counts
        self.pw = pw
        self.wl = wl
        self.NA = NA
        self.n = n
        self.n0 = n0
        self.M = M

        self.Sfloor = Sfloor
        self.alpha = alpha
        self.sigma = sigma
        
        self.pinit = pinit
        self.model = model
        
        self.npix = shape(counts)[0]

        # Pixel coordinate vector
        self.posvec = arange(-(self.npix-1.0)/2.0, (self.npix)/2.0, 1.0) * pw

        # Create instance of dipole PSF based on selected model
        if self.model == 'mortensen':
            self.dd = dipdistr(self.wl, self.n, self.n0, self.M, self.NA)
        elif self.model == 'hinterer':
            self.dd = dipdistr_hinterer(self.wl, self.n, self.n0, self.M, self.NA)
        else:
            raise ValueError(f"Unknown model: {model}. Choose 'mortensen' or 'hinterer'.")
    
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

        # Assign parameters their current values
        mux, muy, b, N, phi, theta, dz = x

        # Convert parameters
        b = b**2
        N = N**2

        # Calculate probabilities for all pixels
        for i in range(npix):
            for j in range(npix):
                if self.model == 'mortensen':
                    pij[j, i] = self.dd.PSF_approx(posvec[i]-mux, posvec[j]-muy, phi, theta, dz)
                else:  # hinterer
                    pij[j, i] = self.dd.PSF_approx(posvec[i]-mux, posvec[j]-muy, phi, theta, dz)
                    
        pij *= (self.pw**2)

        # Subtract noise floor
        effcounts = counts - Sfloor

        # Calculate log-likelihood
        value = 0.0
        for i in range(npix):
            for j in range(npix):
                eta = N*pij[j, i] + b
                
                f0 = alpha * exp(-eta) * eta
                fp0 = f0 * 0.5 * alpha * (eta-2)

                cij = effcounts[j, i]

                # Functions for approximate convolution of the EMCCD readout-noise distribution
                # with the signal distribution from the amplification
                conv0 = 0.5 * (1 + erf(cij/(sqrt(2*sigma**2))))
                conv1 = sigma * exp(-cij**2/(2*sigma**2))/sqrt(2*pi) + cij * conv0
                temp = (f0*conv0 + fp0*conv1 + exp(-eta)*gauss(cij, sigma))
                
                if (cij > 0.0):
                    nij = alpha * cij
                    # Use log-transform for large arguments
                    if eta*nij > 10**5:
                        transform = 0.5*log(alpha*eta/cij) - nij - eta + 2*sqrt(eta*nij) \
                                  - log(2*sqrt(pi)*(eta*nij)**0.25)
                        temp += (exp(transform) - f0 - fp0*cij)
                    else:
                        temp += (sqrt(alpha*eta/cij) * exp(-nij-eta) * i1(2*sqrt(eta*nij)) \
                               - f0 - fp0*cij)

                value += log(temp)

        # The negative log-likelihood is to be minimized
        value *= -1.0
        print(f"{value:7.5f} {mux:5.3f} {muy:5.3f} {b:5.3f} {N:5.3f} {theta:5.3f} {phi:5.3f} {dz:5.3f}")
        return value
