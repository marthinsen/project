#!/usr/bin/env python
"""
   File: a.py
 Author: Eirik Marthinsen
   Date: 2011-03-10 19:49:07 CET

Description: Calculate the angular distribution of scattered light
on a randomly rough surface with pertubation theory.
"""

from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import quad, Inf, quadrature

# Useful constants
c       = 3e8	                                # Speed of light [m/s]
deg2rad	= pi/180	                        # Convert from degrees to radians

# Incoming wave
theta0  = 40 * deg2rad	                        # Angle of incidence [rad]
lambd	= 457.9e-9		                        # Wavelangth [m]
omega   = 2*pi*c/lambd		                # Frequancy of incoming wave [Hz]
K 		= omega/c * sin(theta0)              # Wave number parallell to x1 [1/m]

# Scattered Wave
N       = 20
Theta_s = linspace(-80, 80, num=N)*deg2rad	# Scattering angles [rad]
Q		= omega/c * sin( Theta_s )			# Scattering wave number parallell to x1 [1/m]
I11		= zeros(N)	                        # Array for the 1-1 moments of
I22     = zeros(N)
I31     = zeros(N)

# Define constants
eps = (-7.5 + 0.24j)                        # Permittivity of reflective material [1]
sigma   = 5e-9		                        	# Rms-height of surface [m]
a 		= 100e-9		                        # Transverse correlation distance [m]


def alpha(arg):     # Eq. (3.8)	 [M&M]
    return sqrt(eps * (omega/c)**2 - arg**2)

def alpha0(arg):    # Eq. (3.3)  [M&M]
    return sqrt(0j + (omega/c)**2 - arg**2)

def power(arg):     # Eq. (2.7)  [M&M]
    return sqrt(pi) * a \
            * exp( -a**2 * arg**2 / 4 )

def G0(arg):        # Eq. (3.10) [M&M]
    return 1j * eps / (eps * alpha0(arg) + alpha(arg))

def G02(arg):       # Eq. (3.16) [M&M]
    return (c/omega)**2 * abs(eps)**2 \
            / abs( eps*cos(arg) + sqrt( eps - sin(arg)**2 ) )**2

def Ixx(scAng,Txx): # Eq. (7) [O'D]
  
  return Txx * G02(theta0)*G02(scAng) \
            * 2*omega**3/pi/c**3 * cos(theta0) * cos(scAng)**2

def V1(q,k):      # Eq. (3.20a) [M&M] & (4) [O'D]
    return 1j / eps**2 * (eps - 1) * (eps*q*k - alpha(q)*alpha(k))

def V2x(q,p,k): 
    return 1j * 2*p / (q + k if q + k else 1) \
            * (eps-1) / eps**2 \
            * (alpha(q) + alpha(k)) \
            * (q*k - alpha(q)*alpha(k)) \
           + 2j * (eps-1)**2 / eps**3 \
            * alpha(q) * alpha(p) * alpha(k)

def V2(q,p,k):
    return 1j * (eps - 1) / eps**2 \
            * G(1,q,k) * 2 * p / (q + k if q+k else 1) \
           - (eps - 1) / eps \
            * F(0,q,p) * V1(p,k) 

def G(n,q,k):
    return ( (alpha(q) - alpha0(k))**n \
            * (q*k + alpha(q)*alpha0(k)) \
            * (eps*alpha0(k) - alpha(k)) \
           + (alpha(q) + alpha0(k))**n \
            * (q*k - alpha(q)*alpha0(k)) \
            * (eps*alpha0(k) + alpha(k))) \
           / (2 * alpha0(k) if alpha0(k) else 1)

def F(n,q,p):
    return ( (alpha(q) - alpha0(p))**n \
            * (q*p + alpha(q)*alpha0(p)) \
           - (alpha(q) + alpha0(p))**n \
            * (q*p - alpha(q)*alpha0(p)) ) \
           / (2 * alpha0(p) if alpha0(p) else 1)

def V3(q,p,r,k):
    term1 = 1j * (eps-1) / eps**2 \
             * G(2,q,k) \
             * 6 * r * p * (p*q - r*k) \
             / (q + r if q + r else 1) \
             / (p + k if p + k else 1) \
             / (q**2 - k**2 if q**2 - k**2 else 1)
            # * 6 * p * r \
            # / (q + 2*k if q + 2*k else 1) \
            # / (p + r   if p + r   else 1) \
    term2 = - (eps - 1) / eps \
             * 3 * F(1,q,r) * V1(r,k) \
             * 2 * p / (q + r if q + r else 1)
    term3 = - (eps - 1) / eps \
             * 3 * F(0,q,p) * V2(p,r,k)
    return term1 + term2 + term3        

def V3x(q,p,r,k):
    term1 = 1j * (eps-1) / eps**3 \
            * ( alpha(q) * alpha(k) \
            * ( (eps-2) \
               * (alpha(q)**2 + alpha(k)**2) \
               - 2 * eps * alpha(q) * alpha(k) \
               + 3/2 * (eps-1) * (q**2 + k**2) ) \
            + q * k * ( 2 * alpha(q) * alpha(k) \
               + eps / 2 * (eps + 1) \
                * (alpha(q)**2 + alpha(k)**2) \
               - eps / 2 * (eps - 1) \
                * (q**2 + k**2) ) ) \
            * 6 * r * p * (p*q - r*k) \
            / (q + r if q + r else 1) \
            / (p + k if p + k else 1) \
            / (q**2 - k**2 if q**2 - k**2 else 1)
    term2 = - 6j * (eps - 1)**2 / eps**3 \
            * alpha(r) * alpha(k) * (q*r - alpha(q)**2) \
            * p / (q + r if q + r else 1)
    term3 = - 6j * (eps - 1)**2 / eps**3 \
            * alpha(q) * alpha(p) * (p*k - alpha(k)**2) \
            * r / (p + k if p + k else 1)
    term4 = + 3j * (eps - 1)**2 / eps**3 \
            * ( eps * q * k - alpha(q)*alpha(k) ) \
            * ( r**2 * p / (q + r if q + r else 1) \
              + p**2 * r / (p + k if p + k else 1) )
    term5 = - 6j * (eps - 1 )**3 / eps**4 \
            * alpha(q) * alpha(p) * alpha(r) * alpha(k)
    return term1 + term2 + term3 + term4 + term5

def A2x(q, p, k):
    return 2 * V1(q,p) * G0(p) * V1(p,k)

def A2(q, p, k):
    return V2(q,p,k) \
            + 2 * V1(q,p) * G0(p) * V1(p,k)

def A3x(q, p, r, k):
    term1 = V3(q,p,r,k)
    term2 = 3 * V2(q,p,r) * G0(r) * V1(r,k)
    term3 = 3 * V1(q,p)   * G0(p) * V2(p,r,k)
    term4 = 6 * V1(q,p)   * G0(p) * V1(p,r)   * G0(r) * V1(r,k)
    return term1 + term2 + term3 + term4

def A3(q, p, r, k):
    term1 = V3(q,p,r,k)
    term2 = 3 * V2(q,p,r) * G0(r) * V1(r,k)
    term3 = 3 * V1(q,p)   * G0(p) * A2(p,r,k)
    return term1 + term2 + term3  

def A3xx(q,p,r,k):
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
                       

def integrandT22(p, q, k):
    A2a = A2x(q, p, k)
    A2b = A2x(q, q + k - p, k)
    return A2a \
            * (A2a + A2b).conjugate() \
            * power(q - p) \
            * power(p - k)

def A31(q, p, k):
    term1 = A3xx(q, p+q,   q, k)
    term2 = A3xx(q, p+q, p+k, k)
    term3 = A3xx(q,   k, p+k, k)
    return term1 + term2 + term3

def A311A1x(q,k):
    intgrnd = lambda p: ( (A31(q, p, k) + A31(q, -p, k)) \
                            * V1(q,k).conjugate() ).real \
                            * power(p)
    sing   = omega/c * sqrt(eps/(eps+1)).real
    d      = 1e6
    lims   = [0, sing+k-d, sing+k+d, sing-k-d, sing-k+d, sing-q-d, sing-q+d, sing+q-d, sing+q+d, 1e9,  Inf] 
    lims.sort()
    intgrl = 0 

    for i in xrange(len(lims) - 1):
        print '%3i from %3e to %3e' % (i, lims[i], lims[i+1])
        intgrl += quad(intgrnd, lims[i], lims[i+1])[0]

    return intgrl

def printData():    
    print '*** The Wave ***'
    print '  Wavelength:            %3.1f nm'		% (lambd*1e9)
    print '  Frequency:         %3.3e Hz'			% omega
    print '  Angle in.                 %.0f deg'	% (theta0/deg2rad)
    print '  Wavenum. in:       %3.3e 1/m'			% k
    print '*** The Surface ***'
    print '  rms-height:              %1.1f nm'		% (sigma*1e9)
    print '  Trans. corr. len:        %.0f nm'		% (a*1e9)
    print '  Permittivity:     %1.1f+%.2fi'			% (eps.real,eps.imag)

def iterate():
    for i in xrange(N):
        print '%3i/%i' % (i,N) 
        q       = Q[i]
        theta_s	= Theta_s[i]
    
        T11		= abs(V1(q,K))**2 * sigma**2 * power(q-K)
        I11[i]  = Ixx(theta_s, T11)
    
        T22     = sigma**4 \
                   / (2 * pi) \
                   * ( quad(integrandT22,  -1e10, -0.2e8, args=(q,K))[0] \
                     + quad(integrandT22, -0.2e8, -0.1e8, args=(q,K))[0] \
                     + quad(integrandT22, -0.1e8,  0.1e8, args=(q,K))[0] \
                     + quad(integrandT22,  0.1e8,  0.2e8, args=(q,K))[0] \
                     + quad(integrandT22,  0.2e8,   1e10, args=(q,K))[0] )
        I22[i]  = Ixx(theta_s, T22 / 4)

        T31    = sigma**4 \
                  / (2 * pi) \
                  * A311A1x(q, K) \
                  * power(q - K)
        I31[i] = Ixx(theta_s, - T31 / 3)

def plotIt():
    plot(Theta_s/deg2rad, I11, 'r', \
         Theta_s/deg2rad, I22, 'b', \
         Theta_s/deg2rad, I31, 'g')
    show()
    
def writeToFile():    
    File = open('i.dat','w')
    for i in xrange(N):
        File.write('%e\t%e\t%e\t%e\n' %(Theta_s[i]/deg2rad, \
                                    I11[i], \
                                    I22[i], \
                                    I31[i]))
    File.close()

def foo():
    p = 1e4
    k = 0

    fun1 = lambda q: ((A3(q,p+q,q,k)+A3(q,p+q,p+k,k) \
                      + A3(q,k,p+k,k)) * V1(q,k).conjugate()).real
   
    x1 = linspace(-1e-3, 1e-3,num=1001)

    y1 = array([fun1(i) for i in x1])
   
    plot(x1, y1, 'b')
    show()

def bar():
    color = ['r', 'b', 'g', 'c']
    for k in [0,1,1e5,1e7]:
        q = 0
        I = lambda p: ( (A31(q, p, k) + A31(q, -p, k)) * V1(q,k).conjugate()).real * power(p)
        
        s = omega/c*sqrt(eps/(eps+1))
        x = linspace(0,4e7, num=5000)
        plot(x,map(I,x), color[i])
        axvline(x=s, color='r')
        axvline(x=s+k, color='g')
        axvline(x=s+k+1e6, color='g')
        axvline(x=s+k-1e6, color='g')
        axvline(x=s-k, color='g')
        axvline(x=s-k+1e6, color='g')
        axvline(x=s-k-1e6, color='g')
        axvline(x=s+q, color='c')
        axvline(x=s+q+1e6, color='c')
        axvline(x=s+q-1e6, color='c')
        axvline(x=s-q, color='c')
        axvline(x=s-q-1e6, color='c')
        axvline(x=s-q+1e6, color='c')
    show()

if __name__ == "__main__":
    iterate()
    writeToFile()
    plotIt()
    #foo()
    #bar()
