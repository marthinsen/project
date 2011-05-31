/*
 * =====================================================================================
 *
 *       Filename:  diffscatt.cpp
 *
 *    Description:  Computes the mean diffusive intensity scattered by a rough surface 
 *    				using exact perturbative terms up to fourth order.
 *
 *
 *        Version:  1.0
 *        Created:  31/05/11 12:04:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eirik Marthinsen (eirikma@stud.ntnu.no), 
 *        Company:  NTNU, Trondheim
 *
 * =====================================================================================
 */
#include "diffscatt.h"
#include <gsl/gsl_integration.h>
#include <iostream>
#include <complex>
#include <math.h>
#include <vector>
#include <stdio.h>


// Define the constant c and imaginary number im (i)
#define c 299792458
#define im complex<double>(0,1)

// Incoming wave
double omega  = 2*M_PI*c/457.9e-9;

// Surface constants 
complex<double> eps (-7.5,0.24);
double          sigma = 5e-9;
double          a     = 100e-9;
double 			km    = 0.82 * omega/c;
double 			kp    = 1.29 * omega/c;

// Surface Correlation function
int corr = 1;


double power(double k)
{
	if(corr == 1)
		return sqrt(M_PI)*a*exp(-k*k*a*a/4);
	else
	{
	  	printf(".");
		if( (-kp < k && k< -km) || (km < k && k < kp) )
			return M_PI/(kp+km);
		else
			return 0;
	}
}

complex<double> alpha0(double k)
{
  	complex<double> tmp = pow(omega/c,2) - pow(k,2);
	return sqrt(tmp);
}

complex<double> alpha(double k)
{
	complex<double> tmp = eps*pow(omega/c,2) - pow(k,2);
	return sqrt(tmp);
}

complex<double> G0(double k)
{
	return im*eps/(eps*alpha0(k)+alpha(k));
}

complex<double> A1(double q, double k)
{
	return im*(eps-1.)/pow(eps,2)*(eps*q*k-alpha(q)*alpha(k));
}

complex<double> A2(double q, double p, double k)
{
	return 2.*A1(q,p)*G0(p)*A1(p,k);
}

complex<double> A3(double q, double p, double r, double k)
{
	complex<double> ak, ap, ar, aq;
	ak = alpha(k); ap = alpha(p); ar = alpha(r); aq = alpha(q);
	complex<double> term1, term2, term3, term4, term5, term6, term7;

	term1 = 3.*im * pow(eps-1.,2) / 2. / pow(eps,3)
	         * ( (pow(p,2) + pow(r,2)) * (eps*q*k - aq*ak) 
	                - 2.*(p*k - pow(ak,2))*aq*ap );
	term2 = im * (eps - 1.) / pow(eps,3) *aq*ak 
	         * ( 3./2. * (eps - 1.) * (pow(q,2) + pow(k,2)) 
	                - 2. * eps * aq * ak
	                + (eps - 2.) * (pow(aq,2) + pow(ak,2)) );
	term3 = im * (eps - 1.) / pow(eps,3) * q*k
	         * ( 2.*aq*ak
	                - eps * (eps - 1.) * (pow(q,2) + pow(k,2)) / 2.
	                + eps * (pow(aq,2) + pow(ak,2)) );
	term4 = -3.*im * pow(eps - 1.,2) / pow(eps,3) * ar * ak
	         * (q*r - pow(aq,2));
	term5 = -6.*im * pow(eps - 1.,3) / pow(eps,4) *aq*ap*ar*ak;
	term6 = -3. * pow(eps - 1.,2) / pow(eps,4)
	         * (eps*r*k - ar*ak) * G0(r)
	         * ( 2. * (eps - 1.) / eps * aq * ap * ar
	                + (aq + ar) * (q*r - aq*ar) );
	term7 = -3. * pow(eps - 1.,2) / pow(eps,4)
	         * (eps*q*p - aq*ap) * G0(p)
	         * ( 2. * (eps - 1.) / eps * ap*ar*ak
	                + (ap + ak) * (p*k - ap*ak)
	                + 2.*im * (eps - 1.) / pow(eps,2)
	                * (eps*p*r - ap*ar) * G0(r)
	                * (eps*r*k - ar*ak) );
	return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

double integrandA311Re(double p, void * params)
{
    double *param = (double *)params;
	double q = *(double *) param++;
	double k = *(double *) param;
	return  real(A3(q,p+q,q,k) + A3(q,p+q,k+p,k) + A3(q,k,p+k,k)) * power(p);
}

double integrandA311Im(double p, void * params)
{
    double *param = (double *)params;
	double q = *(double *) param++;
	double k = *(double *) param;
	return  imag(A3(q,p+q,q,k) + A3(q,p+q,k+p,k) + A3(q,k,p+k,k)) * power(p);
}

complex<double> A311(double q, double k)
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	double params[] = {q,k};
	double IntRe = 0; double IntIm = 0;
    double lim[] = {-1e11, -1e8, 1e7, 1e6, 0, 1e6, 1e7, 1e8, 1e11};
	int  lims = sizeof(lim)/sizeof(double);
	gsl_function funcRe;
	funcRe.function = &integrandA311Re;
	funcRe.params   = &params;
	gsl_function funcIm;
	funcIm.function = &integrandA311Im;
	funcIm.params   = &params;
	for(int j=1; j<lims; ++j) 
	{
		gsl_integration_qags (&funcRe, lim[j-1], lim[j], 0, 1e-7, 1000, w, &result, &error);
		IntRe += result;
		gsl_integration_qags (&funcIm, lim[j-1], lim[j], 0, 1e-7, 1000, w, &result, &error);
		IntIm += result;
	}
	return IntRe + im*IntIm;
}

double integrandT22(double p, void * params)
{
    double *param = (double *)params;
	double q = *(double *) param++;
	double k = *(double *) param;
	return real(A2(q,p,k)*conj(A2(q,p,k) + A2(q,q+k+p,k))*power(q-p)*power(p-k));
}

double Ixx(double q, double k, double Txx)
{
  	double v0 = asin(k*c/omega);
	double vs = asin(q*c/omega);
	return  2/M_PI *pow(omega/c,3) *pow(cos(vs),2) *cos(v0) *pow(abs(G0(k)),2) *pow(abs(G0(q)),2) *Txx;
}

void iterate11 (vector<double>* I11, int n, double theta, double min, double max)
{
	double k = omega/c*sin((theta)*M_PI/180);
	for(int i=0; i<n; ++i)
	{
		double q = omega/c*sin((-min+i*(max-min)/(n-1))*M_PI/180);
		double T11 = pow(sigma,2)*abs(pow(A1(q,k),2))*power(q-k);
		I11->push_back( Ixx(q, k, T11) );
	}
	return;
}

void iterate22 (vector<double>* I22, int n, double theta, double min, double max)
{
	double k = omega/c*sin((theta)*M_PI/180);
	for(int i=0; i<n; ++i)
	{
		printf("%4.1f %%\n", 1.*i*100/n); 
		double q = omega/c*sin((-min+i*(max-min)/(n-1))*M_PI/180);
		double lim[] = {-1e10, -2e7, -1e7, 1e7, 2e7, 1e10};
		int  lims = sizeof(lim)/sizeof(double);
		double Int = 0;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	    double result, error;
		double params[] = {q,k};
		gsl_function FUNC;
	    FUNC.function = &integrandT22;
		FUNC.params   = &params;
	    for(int j=1; j<lims; ++j) 
		{
			gsl_integration_qags (&FUNC, lim[j-1], lim[j], 0, 1e-7, 1000, w, &result, &error);
			Int += result;
		}
		double T22 = pow(sigma,4) /2/M_PI * Int;
		I22->push_back( Ixx(q, k, T22/4) );
	}
}

void iterate31 (vector<double>* I31, int n, double theta, double min, double max)
{
	double k = omega/c*sin((theta)*M_PI/180);
	for(int i=0; i<n; ++i)
	{
		printf("%4.1f %%\n", 1.*i*100/n); 
		double q = omega/c*sin((-min+i*(max-min)/(n-1))*M_PI/180);
		double T31 = pow(sigma,4)/2./M_PI*real(conj(A1(q,k))*A311(q,k))*power(q-k);
		I31->push_back( Ixx(q, k, T31/3));
	}
	return;
}
