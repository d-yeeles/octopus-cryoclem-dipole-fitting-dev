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


import numpy as np


class ChirpZTransform:
    """
    2D chirp-z transform

    See the following publication:
    Raoqiong Tong and Robert W. Cox.
    "Rotation of NMR images using the 2D chirp-z transform"
    Magnetic Resonance in Medicine 41.2 (1999): 253-256.
    """

    def __init__(self):

        # Get dimensions of image and k-space
        nPixels = 19 # !!! remove this hardcoded value later !!!
        oversampling = 3 # !!! remove this hardcoded value later !!!
        nDiscretizationBFP = 129 # !!! remove this hardcoded value later !!!
        pixelSize = 51.2 # !!! remove this hardcoded value later !!
        nPixelsImage = nPixels * oversampling  # size of field in pixels
        nPixelsBFP = nDiscretizationBFP  # size of input field (k-space)
        self.nPixelsPadded = nPixelsImage + nPixelsBFP - 1  # padded size required for convolution

        NA = 1.49
        wl = 580

        unitObjectSpace = (pixelSize/1e9) / oversampling
        # Largest spatial frequency passed by objective lens
        maxSpatialFrequency = NA / (wl/1e9)
        maxAngularFrequency = 2*pi * maxSpatialFrequency
        # Unit in pupil space (k-space)
        unitKSpace = 2 * maxAngularFrequency / nDiscretizationBFP # full range covers 2 times the maximum frequency




        # Create grid
        x = np.arange(np.floor(-self.nPixelsPadded / 2 + 0.5), np.floor(self.nPixelsPadded / 2 - 0.5) + 1)
        y = np.arange(np.floor(-self.nPixelsPadded / 2 + 0.5), np.floor(self.nPixelsPadded / 2 - 0.5) + 1)
        X, Y = np.meshgrid(y, x)

        # Create kernel (quadratic convolution phase kernel)
        alpha = unitObjectSpace * unitKSpace / (2 * np.pi)
        self.kernel = np.exp(-1j * alpha * np.pi * (X ** 2 + Y ** 2))

        # Fourier transform of kernel
        self.fourierKernel = np.fft.fft2(np.fft.ifftshift(self.kernel))

    def apply(self, E_in):

        nPixels = 19 # !!! remove this hardcoded value later !!!
        oversampling = 3 # !!! remove this hardcoded value later !!!

        # Complex conjugate of kernel
        conjKernel = np.conj(self.kernel)

        # Electric field times complex conjugate of kernel
        f = self.embed_array_2D(E_in, self.nPixelsPadded, 0) * conjKernel

        # Fourier transform of kernel
        F = np.fft.fft2(np.fft.ifftshift(f))

        # Convolving f with kernel (inverse fft of multiplication)
        convolution = np.fft.fftshift(np.fft.ifft2(F * self.fourierKernel))

        # Final multiplication by factor
        E_out = self.nPixelsPadded * conjKernel * convolution

        # Crop to image size
        E_out = self.crop_array_2D(E_out, nPixels * oversampling)

        return E_out

    @staticmethod
    def embed_array_2D(array, new_size, pad_value):
        """Embeds a smaller array into a larger zero-padded array."""
        padded_array = np.full(new_size, pad_value, dtype=array.dtype)
        insert_slices = tuple(slice((ns - s) // 2, (ns - s) // 2 + s) for s, ns in zip(array.shape, new_size))
        padded_array[insert_slices] = array
        return padded_array

    @staticmethod
    def crop_array_2D(array, new_size):
        """Crops a larger array to a smaller size."""
        crop_slices = tuple(slice((s - ns) // 2, (s - ns) // 2 + ns) for s, ns in zip(array.shape, new_size))
        return array[crop_slices]


class Mask:
    def __init__(self, par=None):
        self.values = None
        self.nGrid = 129
        self.spacing = 1  # spacing between two neighboring entries
        self.mode = 'FFT'  # 'exact' or 'FFT'
        self.radius = None

        if par:
            self.set_input_parameters(par)

        self.calculate_mask()

    def set_input_parameters(self, par):
        for key, value in par.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def calculate_mask(self):
        N = [self.nGrid, self.nGrid]
        s = [self.spacing, self.spacing]

        if self.mode == 'exact':
            x = np.linspace(-(N[0] - 1) / 2, (N[0] - 1) / 2, N[0]) * s[0]
            y = np.linspace(-(N[1] - 1) / 2, (N[1] - 1) / 2, N[1]) * s[1]
        elif self.mode == 'FFT':
            x = np.arange(np.floor(-N[0] / 2 + 0.5), np.floor(N[0] / 2 - 0.5) + 1) * s[0]
            y = np.arange(np.floor(-N[1] / 2 + 0.5), np.floor(N[1] / 2 - 0.5) + 1) * s[1]

        X, Y = np.meshgrid(y, x)
        self.radius = np.sqrt(X ** 2 + Y ** 2)
        self.values = (self.radius ** 2 <= (min(np.array(N) / 2 * np.array(s))) ** 2)

    def get_polar_coordinates(self):
        x = np.linspace(-1, 1, self.nGrid)
        X, Y = np.meshgrid(x, x)
        angle, normalized_radius = np.arctan2(Y, X), np.sqrt(X ** 2 + Y ** 2)
        angle = np.fliplr(angle + np.pi)

        return normalized_radius, angle

    def plot(self):
        plt.imshow(self.values, cmap='gray', origin='lower')
        plt.axis('equal')
        plt.axis('tight')
        plt.show()

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

    def Intensity_approx(self,x,y,phi,theta,deltaz):
        # # Transform alpha to define optical axis in direction from camera towards objective
        # alpha+=pi
        # self.deltaz=deltaz
        # r=sqrt(x**2+y**2)
        # if (x<0.0): phip=pi-arctan(-y/x)
        # else: phip=arctan(y/x)
        # rho=self.M*r

        # Calculate electric field
        RI = [self.np, self.n, self.n0]
        h_intermediate = 1 # not what this is
        pos = [x, y, deltaz]
        z = pos[2]
        diameter_of_pupil = 1e-3
        focal_length = self.NA*diameter_of_pupil
        other_mu = 1e-12

        # Pre-Calculations
        if len(RI) == 1:
            RI = [RI[0], RI[0], RI[0]]

        # Coordinates in the objective pupil
        pupil_mask = Mask()
        _, mask_angle = pupil_mask.get_polar_coordinates()
        PHI3 = np.fliplr(mask_angle) - np.pi

        # Wavenumbers (magnitude of k-vectors) in different media
        k0 = 2 * np.pi / (self.wl/1e9)
        k1 = k0 * RI[0]
        k2 = k0 * RI[1]  # Wavenumber in media 2 (intermediate layer)
        k3 = k0 * RI[2]  # Wavenumber in media 3 (immersion medium)

        # Angles in different media
        Kr = pupil_mask.radius
        THETA1 = np.arccos(np.sqrt(1 - (RI[2] / RI[0] * Kr / k3) ** 2)) * pupil_mask.values
        THETA2 = np.arccos(np.sqrt(1 - (RI[2] / RI[1] * Kr / k3) ** 2)) * pupil_mask.values
        THETA3 = np.arcsin(Kr / k3) * pupil_mask.values

        # Cosines of angles
        CTHETA1 = np.cos(THETA1)
        CTHETA2 = np.cos(THETA2)
        CTHETA3 = np.cos(THETA3)

        # Fresnel coefficients
        tp12 = (2 * RI[0] * CTHETA1) / (RI[0] * CTHETA2 + RI[1] * CTHETA1)
        tp23 = (2 * RI[1] * CTHETA2) / (RI[1] * CTHETA3 + RI[2] * CTHETA2)
        ts12 = (2 * RI[0] * CTHETA1) / (RI[0] * CTHETA1 + RI[1] * CTHETA2)
        ts23 = (2 * RI[1] * CTHETA2) / (RI[1] * CTHETA2 + RI[2] * CTHETA3)
        rp12 = (RI[1] * CTHETA1 - RI[0] * CTHETA2) / (RI[0] * CTHETA2 + RI[1] * CTHETA1)
        rp23 = (RI[2] * CTHETA2 - RI[1] * CTHETA3) / (RI[1] * CTHETA3 + RI[2] * CTHETA2)
        rs12 = (RI[0] * CTHETA1 - RI[1] * CTHETA2) / (RI[0] * CTHETA1 + RI[1] * CTHETA2)
        rs23 = (RI[1] * CTHETA2 - RI[2] * CTHETA3) / (RI[1] * CTHETA2 + RI[2] * CTHETA3)

        # Fresnel coefficients for three-layer system
        tp = (tp12 * tp23 * np.exp(1j * k2 * h_intermediate * CTHETA2)) / (
                    1 + rp12 * rp23 * np.exp(2j * k2 * h_intermediate * CTHETA2))
        ts = (ts12 * ts23 * np.exp(1j * k2 * h_intermediate * CTHETA2)) / (
                    1 + rs12 * rs23 * np.exp(2j * k2 * h_intermediate * CTHETA2))

        # Dipole projections
        mu_p = other_mu * np.sin(theta) * np.cos(phi - PHI3)
        mu_s = other_mu * np.sin(theta) * np.sin(phi - PHI3)
        mu_z = other_mu * np.cos(theta)

        # Prefactor C
        C = ((k3 ** 2 * np.exp(1j * k3 * focal_length) * CTHETA3) / (focal_length * RI[0])) * np.exp(
            -1j * k3 * h_intermediate * CTHETA3) * np.exp(1j * k1 * CTHETA1 * z)

        # Electric field components
        E3p = C * tp * CTHETA3 * (mu_p / RI[2] + mu_z * np.sin(THETA3) / CTHETA1)
        E3s = C * ts * (mu_s / (RI[2] / CTHETA1))
        E3z = C * tp * np.sin(THETA3) * (mu_p / RI[2] + mu_z * np.sin(THETA3) / CTHETA1)
        # Apodization factor
        apodization_factor = (1 / np.sqrt(CTHETA3)) * pupil_mask.values
        E_BFP_p = (E3p * CTHETA3 + E3z * np.sin(THETA3)) * apodization_factor
        E_BFP_s = E3s * apodization_factor

        # Coordinate transformation
        E_BFP_x = np.cos(PHI3) * E_BFP_p - np.sin(PHI3) * E_BFP_s
        E_BFP_y = np.sin(PHI3) * E_BFP_p + np.cos(PHI3) * E_BFP_s

        # get intensity
        # cz_transform = ChirpZTransform()
        # Ix = abs(cz_transform.apply(E_BFP_x))**2 # intensity = | E_imagePlane |²
        # Iy = abs(cz_transform.apply(E_BFP_y))**2 # intensity = | E_imagePlane |²
        Ix = abs(E_BFP_x)**2 # intensity = | E_imagePlane |²
        Iy = abs(E_BFP_y)**2 # intensity = | E_imagePlane |²

        psf = Ix + Iy

        return psf

    def PSF_approx(self,x,y,alpha,beta,deltaz):
        # Transform alpha to define optical axis in direction from camera towards objective
        alpha+=pi
        self.deltaz=deltaz
        r=sqrt(x**2+y**2)
        if (x<0.0): phip=pi-arctan(-y/x)
        else: phip=arctan(y/x)
        rho=self.M*r
        value=self.M**2*self.Intensity_approx(x,y,alpha,beta,deltaz)
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

    def Value(self,x):

        """ Calculation of the negative log-likelihood

            Input:
            ------------------------------------------------------------
            x (array)       : Array of current parameter values

            Output:
            ------------------------------------------------------------
            value (float)   : The negative log-likelihood
        """

        counts=self.counts
        npix=self.npix
        posvec=self.posvec

        Sfloor=self.Sfloor
        alpha=self.alpha
        sigma=self.sigma

        pij=zeros((npix,npix))

        params=zeros(7)

        # Assign parameters their current values
        mux,muy,b,N,phi,theta,dz=x

        # Convert parameters
        b=b**2
        N=N**2

        # Calculate probabilities for all pixels
        for i in range(int(npix)):
            for j in range(int(npix)):
                pij[j,i]=self.dd.PSF_approx(posvec[i]-mux,posvec[j]-muy,phi,theta,dz)
        pij*=(self.pw**2)

        # Subtract noise floor
        effcounts=counts-Sfloor

        # Calculate log-likelihood
        value=0.0
        for i in range(int(npix)):
            for j in range(int(npix)):
                eta=N*pij[j,i]+b

                f0=alpha*exp(-eta)*eta
                fp0=f0*0.5*alpha*(eta-2)

                cij=effcounts[j,i]

                # Functions for approximate convolution of the EMCCD readout-noise distribution
                # with the signal distribution from the amplification
                conv0=0.5*(1+erf(cij/(sqrt(2*sigma**2))))
                conv1=sigma*exp(-cij**2/(2*sigma**2))/sqrt(2*pi)+cij*conv0
                temp=(f0*conv0+fp0*conv1+exp(-eta)*gauss(cij,sigma))

                if (cij>0.0):
                    nij=alpha*cij
                    # Use log-transform for large arguments
                    if eta*nij>10**5:
                        transform=0.5*log(alpha*eta/cij)-nij-eta+2*sqrt(eta*nij)\
                                   -log(2*sqrt(pi)*(eta*nij)**0.25)
                        temp+=(exp(transform)-f0-fp0*cij)
                    else:
                        temp+=(sqrt(alpha*eta/cij)*exp(-nij-eta)*i1(2*sqrt(eta*nij))\
                            -f0-fp0*cij)

                # dave feb 2025 - temp sometimes hits 0 due to rounding errors
                if temp <= 0:
                    temp = 1e-6

                value+=log(temp)

        # The negative log-likelihood is to be minimized
        value*=-1.0
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

        # Transformation of initial values
        pinit=zeros(7)
        pinit[0:2]=self.initvals[0:2] # x, y
        pinit[2]=sqrt(self.initvals[2]) # b
        pinit[3]=sqrt(self.initvals[3]) # N
        pinit[4]=self.initvals[4] # phi
        pinit[5]=self.initvals[5] # theta
        pinit[6]=self.initvals[6] # dz

        # Create instance of LogLikelihood object
        ll=LogLikelihood(counts,self.a,self.wl,self.NA,self.n,self.n0,self.M,pinit,\
                         self.Sfloor,self.alpha,self.sigma)

        # Perform maximization of the log-likelihood using Powell's method
        start_powell = time.time()
        res=fmin_powell(ll.Value,pinit,ftol=0.000000001,maxiter=15,maxfun=10)
        end_powell = time.time()
        elapsed_time_powell = end_powell - start_powell
        print(f"Time: {elapsed_time_powell:.4f} seconds to maximization log-likelihood")

        est=res#[0]

        # # Perform maximization of the log-likelihood using Powell's method
        # start_powell = time.time()
        # # res=fmin_powell(ll.Value,pinit,ftol=0.000000001,maxiter=50,maxfun=500)
        # bounds = [(-1000, 1000), (-1000, 1000), (None, None), (None, None), (0, 2*np.pi), (0, np.pi/2), (None, None)]
        # options = {
        #     'maxiter': 100,  # Maximum number of iterations
        # }
        # res=minimize(ll.Value,pinit,method='trust-constr',bounds=bounds, options=options)
        # end_powell = time.time()
        # elapsed_time_powell = end_powell - start_powell
        # print(f"Time: {elapsed_time_powell:.4f} seconds to maximization log-likelihood")
        # est=res.x#[0]

        # warnflag=res[5]

        L=ll.Value(est)

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

        # print("\nx coordinate [nm] = ", around(est[0],3))
        # print("y coordinate [nm] = ", around(est[1],3))
        # print("azimuthal angle [rad] = ", around(est[4],2))
        # print("polar angle [rad] = ", around(est[5],2))

        x_est = est[0]
        y_est = est[1]
        phi_est = est[4] % (2 * np.pi)
        theta_est = est[5] % (2 * np.pi)
        # cov_mat = covarmatrix

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


#
# if __name__=='__main__':
#
#     close('all')
#
#     # Experimental parameters
#     wl=580.0
#     n=1.52
#     n0=1.33
#     NA=1.49
#     M=250.0
#
#     deltax=32.0
#
#     # PSF parameters
#     N=10000.0
#     b=5.0
#
#     mu=0.1
#     nu=0.1
#     phi=2*pi/3.0
#     theta=0.5
#
#     # Distance from design focal plane
#     deltaz=-30.0
#
#     # EMCCD parameters
#     inversegain=0.0089
#     sigmanoise=12.6
#     Sfloor=300.0
#     gain=1.0/inversegain
#
#     # Pixel array parameters
#     npix=12
#     deltapix=npix/2
#
#     # Load data
#     datamatrix=numpy.loadtxt('data.txt')
#
#     # Initial guess
#     initvals = array([mu,nu,b,N,phi,theta,deltaz])
#     initpix=(deltapix,deltapix)
#
#     # Create instance of MLEwT
#     track=MLEwT(wl,deltax,M,NA,n,np,n0,initvals,initpix,deltapix,Sfloor,inversegain,sigmanoise)
#     # Perform estimation
#     track.Estimate(datamatrix)
#
#     show()
