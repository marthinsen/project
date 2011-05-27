#!/usr/bin/env python
"""
   File: a.py
 Author: Eirik Marthinsen
   Date: 2011-03-10 19:49:07 CET

Description: Calculate the angular distribution of scattered light
on a randomly rough surface with pertubation theory.
"""

from numpy import *
from matplotlib .pyplot import *
from scipy.integrate import quad, Inf

# Physical constant
c       = 3e8	                                # Speed of light [m/s]


# Incoming wave
theta0  = 20 * pi/180	                        # Angle of incidence [rad]
lambd	= 457.9e-9		                        # Wavelangth [m]
omega   = 2 * pi * c / lambd                    # Frequancy of incoming wave [Hz]


# Scattered Wave
N       = 19 
Theta_s = linspace(-89, 89, num=N) * pi/180   	# Scattering angles [rad]
Q		= omega / c * sin( Theta_s )			# Scattering wave number parallell to x1 [1/m]
I11		= zeros(N)  	                        # Array for the 1-1 moments of
I22     = zeros(N)
I31     = zeros(N)

# Define constants
eps = (-7.5 + 0.24j)                            # Permittivity of reflective material [1]
sigma   = 5e-9		                        	# Rms-height of surface [m]
a 		= 100e-9		                        # Transverse correlation distance [m]

def power(k):
    return sqrt(pi) * a * exp(- k**2 * a**2 / 4)

def alpha0(k):
    return sqrt(0j + omega**2 / c**2 - k**2)

def alpha(k):
    return sqrt(eps * omega**2 / c**2 - k**2)

def G0(k):
    return 1j * eps \
            / (eps * alpha0(k) + alpha(k))

def G02(theta):
    return c**2 / omega**2 * abs(eps)**2 \
            / abs(eps*cos(theta) \
              + sqrt(eps - sin(theta)**2))**2

def F(n, q, p):
    term1 = (alpha(q) - alpha0(p))**n \
            * (q*p + alpha(q)*alpha0(p))
    term2 = (alpha(q) + alpha0(p))**n \
            * (q*p - alpha(q)*alpha0(p))
    denom = 2 * alpha0(p)
    return (term1 - term2) / (denom if denom else 1)

def G(n, q, k):
    term1 = (alpha(q) - alpha0(k))**n \
            * (q*k + alpha(q)*alpha0(k)) \
            * (eps*alpha0(k) - alpha(k))
    term2 = (alpha(q) + alpha0(k))**n \
            * (q*k - alpha(q)*alpha0(k)) \
            * (eps*alpha0(k) + alpha(k))
    denom = 2 * alpha0(k)
    return (term1 + term2) / (denom if denom else 1)

def V1(q, k):
    return 1j * (eps - 1) / eps**2 \
            * (eps*q*k - alpha(q)*alpha(k))

def V2(q, p, k):
    term1 = 1j * (eps-1) / eps**2 * G(1, q, k) \
            * 2 * p / (q + k if q + k else 1)
    term2 = - 2 * (eps - 1) / eps \
            * alpha(q) * alpha(p) *  V1(p, k)
    return term1 + term2
        
def V3(q, p, r, k):
    term1 = 1j * (eps - 1) / eps**2 * G(2,q,k) \
            * 6 * r * p * (p*q - r*k) \
            / (q**2 - k**2 if q**2 - k**2  else 1) \
            / (p + k if p + k else 1) \
            / (q + r if q + r else 1)
    term2 = - (eps - 1) / eps \
            * 6 * (alpha(q)**2 - q*r) * V1(r, k) \
            * p / (q + r if q + r else 1)
    term3 = - (eps - 1) / eps \
            * 3 * alpha(q) * alpha(p) * V2(p,r,k)
    return term1 + term2 + term3

def A2(q, p, k):
    return V2(q, p, k) \
            + 2 * V1(q,p) * G0(p) * V1(p,k)

def A2x(q, p, k):
    return 2 * V1(q,p) * G0(p) * V1(p,k)
        
def A3(q, p, r, k):
    term1 = V3(q, p, r, k)
    term2 = 3 * V2(q, p, r) * G0(r) * V1(r, k)
    term3 = 3 * V1(q, p)    * G0(p) * A2(p, r, k)
    return term1 + term2 + term3

def A3x(q,p,r,k):
    ak = alpha(k); ap = alpha(p); ar = alpha(r); aq = alpha(q)

    term1 = 3j * (eps - 1)**2 / 2 / eps**3 \
            * ( (p**2 + r**2) * (eps*q*k - aq*ak) \
               - 2 * (p*k - ak**2) * (aq*ap) )
    term2 = 1j * (eps - 1) / eps**3 * aq * ak \
            * ( 3. / 2 * (eps - 1) * (q**2 + k**2) \
               - 2 * eps * aq * ak \
               + (eps - 2) * (aq**2 + ak**2) )
    term3 = 1j * (eps - 1) / eps**3 * q * k \
            * ( 2 * aq * ak \
               - eps * (eps - 1 ) * (q**2 + k**2) / 2. \
               + eps * (aq**2 + ak**2) )
    term4 = -3j * (eps - 1)**2 / eps**3 * ar * ak \
            * (q*r -aq**2)
    term5 = -6j * (eps - 1)**3 / eps**4 *aq*ap*ar*ak
    term6 = -3 * (eps - 1)**2 / eps**4 \
            * (eps*r*k - ar*ak) * G0(r) \
            * ( 2 * (eps - 1) / eps * aq * ap * ar \
               + (aq + ar) * (q*r - aq*ar) )
    term7 = -3 * (eps - 1)**2 / eps**4 \
            * (eps*q*p - aq*ap) * G0(p) \
            * ( 2 * (eps - 1) / eps * ap * ar * ak \
               + (ap + ak) * (p*k + ap*ak) \
               + 2j * (eps - 1) / eps**2 \
               * (eps*p*r - ap*ar) * G0(r) \
               * (eps*r*k - ar*ak) )
    return term1 + term2 + term3 + term4 + term5 + term6 + term7
 
def A31(q, p, k):
    term1 = A3x(q, q+p,   q, k)
    term2 = A3x(q, q+p, p+k, k)
    term3 = A3x(q,   k, p+k, k) 
    return term1 + term2 + term3

def A311(q,k):
    fRe = lambda p: (A31(q, p, k) + A31(q, -p, k)).real * power(p)
    fIm = lambda p: (A31(q, p, k) + A31(q, -p, k)).imag * power(p)
    s = (omega/c*sqrt(eps/(eps+1))).real
    d = 5e5
    lim = [0, s-k-d, s-k+d, s-q-d, s-q+d, s+k-d, s+k+d, s+q-d, s+q+d, 5e7, 8e7, 1e8, 1e10, Inf]
    lim.sort()
    iRe = iIm = 0
    for i in xrange(len(lim) - 1):
        print '%2.4e - %2.4e' % (lim[i],lim[i+1])
        iRe += quad(fRe, lim[i], lim[i+1])[0]
        iIm += quad(fIm, lim[i], lim[i+1])[0]
    return iRe + 1j*iIm

def Ixx(thetaI, thetaS, Txx):
    return 2 / pi * omega**3 / c**3 \
            * cos(thetaS)**2 * cos(thetaI) \
            * G02(thetaS) * Txx * G02(thetaI)

def iterate11():
    print 'Iterating over the 11 term'
    k = omega / c * sin(theta0)
    for i in xrange(N):
        q      = Q[i]
        
        T11    = sigma**2 \
                  * abs(V1(q,k))**2 \
                  * power(q-k)

        I11[i] = Ixx(theta0, Theta_s[i], T11)

def iterate22():
    print 'Iterating over the 22 term'
    k = omega / c * sin(theta0)
    for i in xrange(N):
        print '%3.1f %%' % (1.*i/N*100)
        q      = Q[i]
        
        fRe = lambda p: ( A2x(q, p, k) \
                        * (A2x(q, p, k) + A2x(q, q + k - p, k)).conjugate() ) \
                      * power(q - p) * power(p - k)
                        
        lim    = [-1e10, -2e7, -1e7, 1e7, 2e7, 1e10] 
        
        Re = sum([quad(fRe, lim[j], lim[j+1])[0] for j in xrange(len(lim)-1)])

        T22    = sigma**4 / (2*pi) * (Re)\
    
        I22[i] = Ixx(theta0, Theta_s[i], T22 / 4)

def iterate31():
    print 'Iterating over the 31 term'
    k = omega / c * sin(theta0)
    for i in xrange(N):
        q      = Q[i]
        T31    = sigma**4 / (2*pi) \
                  * A311(q, k) \
                  * V1(q, k).conjugate() \
                  * power(q - k)
        I31[i] = Ixx(theta0, Theta_s[i], T31.real / 3)
        print '%3.1f %%' % (1.*(i+1)/N*100)

def plotIt():
    plot(Theta_s*180/pi, I11, \
         Theta_s*180/pi, I22, \
         Theta_s*180/pi, I31)
    show()
    
def writeIt():
    File = open('out3.dat','w')
    for i in xrange(N):
        File.write('%e\t%e\t%e\t%e\n' % (Theta_s[i] / pi * 180, \
                                 I11[i], \
                                 I22[i], \
                                 I31[i]))
    File.close()
    
if __name__ == "__main__":
    #q = 0; p = 0 
    #fr = lambda k: (A31(q, p, k) + A31(q, -p, k)).real * power(p)
    #fi = lambda k: (A31(q, p, k) + A31(q, -p, k)).imag * power(p)
    #
    #x = linspace(-1e-10,1e-10,201)
    #y = [fr(x[i]) for i in xrange(len(x))]
    #z = [fi(x[i]) for i in xrange(len(x))]
    #
    #plot(x,y,x,z)
    #show()
    #
    #print fr(0)
    #print fi(0)
    
    #iterate11()
    iterate22()
    #iterate31()
    
    writeIt()
    plotIt()

    print I31
