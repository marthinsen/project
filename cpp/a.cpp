#include <math.h>
#include <gsl/gsl_integration.h>
#include <vector>
#include <complex>
#include <iostream>
#include <stdio.h>
using namespace std;

// Physical constants
const double          c = 299792458;
const complex<double> im = sqrt( complex<float>(-1) );

// Surface
const complex<double> eps (-7.5,0.24);
const double          sigma = 5e-9;
const double          a     = 100e-9;

// Incoming wave
      double theta0 = 40*M_PI/180;
const double lambda = 457.9e-9;
const double omega  = 2*M_PI*c/lambda;

// Scattered wave
int N = 201;

double power(double k);
complex<double> alpha0(double k);
complex<double> alpha(double k);
complex<double> G0(double k);
complex<double> V1(double q, double k);
complex<double> V2(double q, double p,double k);
complex<double> A2x(double q, double p, double k);
complex<double> A3(double q, double p, double r, double k);
double intA311Re(double p, void * params);
double intA311Im(double p, void * params);
complex<double> A311(double q, double k);
double intT22(double p, void * params);
double Ixx(double q, double k, double Txx);
void readInput(void);

int main(void)
{
  	readInput();

  	double k = omega/c*sin(theta0);
	cout << "Iterating over the 11 term" << endl;
	vector<double> I11; 
	for(int i=0; i<N; ++i)
	{
		double q = omega/c*sin((-89+i*178./(N-1))*M_PI/180);
		double T11 = pow(sigma,2)*abs(pow(V1(q,k),2))*power(q-k);
		I11.push_back( Ixx(q, k, T11) );
	}
		
	cout << "Iterating over the 22 term" << endl;
	vector<double> I22; 
	for(int i=0; i<N; ++i)
	{
		printf("%4.1f %%\n", 1.*i*100/N); 
		double q = omega/c*sin((-89+i*178./(N-1))*M_PI/180);
		double lim[] = {-1e10, -2e7, -1e7, 1e7, 2e7, 1e10};
		int  lims = sizeof(lim)/sizeof(double);
		double Int = 0;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	    double result, error;
		double params[] = {q,k};
		gsl_function FUNC;
	    FUNC.function = &intT22;
		FUNC.params   = &params;
	    for(int j=1; j<lims; ++j) 
		{
			gsl_integration_qags (&FUNC, lim[j-1], lim[j], 0, 1e-7, 1000, w, &result, &error);
			Int += result;
		}
		double T22 = pow(sigma,4) /2/M_PI * Int;
		I22.push_back( Ixx(q, k, T22/4) );
	}

	cout << "Iterating over the 31 term" << endl;
	vector<double> I31; 
	for(int i=0; i<N; ++i)
	{
		printf("%4.1f %%\n", 1.*i*100/N); 
		double q = omega/c*sin((-89+i*178./(N-1))*M_PI/180);
		double T31 = pow(sigma,4)/2./M_PI*real(conj(V1(q,k))*A311(q,k))*power(q-k);
		I31.push_back( Ixx(q, k, T31/3));
	}
	  
	FILE * File;
	File = fopen ("out.dat","w");
	for(int i=0; i<N; ++i)
	{
	  	double theta = -89+i*178./(N-1);
		fprintf(File, "%e\t%e\t%e\t%e\n", theta, I11[i], I22[i], I31[i]);
	}
	fclose (File);
	return 0;
}

double power(double k)
{
	return sqrt(M_PI)*a*exp(-k*k*a*a/4);
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

complex<double> V1(double q, double k)
{
	return im*(eps-1.)/pow(eps,2)*(eps*q*k-alpha(q)*alpha(k));
}

complex<double> V2(double q, double p,double k)
{
  	complex<double> tmp = q+k ? 2*p/(q+k) : 1;
	return im*(eps-1.)/pow(eps,2)*(alpha(q)+alpha(k)) + im*(eps-1.)/pow(eps,3)*alpha(q)*alpha(p)*alpha(k);
}

complex<double> A2x(double q, double p, double k)
{
	return 2.*V1(q,p)*G0(p)*V1(p,k);
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

double intA311Re(double p, void * params)
{
    double *param = (double *)params;
	double q = *(double *) param++;
	double k = *(double *) param;
	return  real(A3(q,p+q,q,k) + A3(q,p+q,k+p,k) + A3(q,k,p+k,k)) * power(p);
}

double intA311Im(double p, void * params)
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
	funcRe.function = &intA311Re;
	funcRe.params   = &params;
	gsl_function funcIm;
	funcIm.function = &intA311Im;
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

double intT22(double p, void * params)
{
    double *param = (double *)params;
	double q = *(double *) param++;
	double k = *(double *) param;
	return real(A2x(q,p,k)*conj(A2x(q,p,k) + A2x(q,q+k+p,k))*power(q-p)*power(p-k));
}

double Ixx(double q, double k, double Txx)
{
  	double v0 = asin(k*c/omega);
	double vs = asin(q*c/omega);
	return  2/M_PI *pow(omega/c,3) *pow(cos(vs),2) *cos(v0) *pow(abs(G0(k)),2) *pow(abs(G0(q)),2) *Txx;
}

void readInput(void)
{
	FILE *input = fopen("par.in", "r");
	if ( input == NULL )
	{
		printf("Unable to read inputfile \'par.in\'");
		return;
	}

	char line[80];
	if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%d", &N) )
		printf("Problem reading frequence.\n");
	else
		printf("Will evaluate %d points.\n", N);

	double v;
	if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &v) )
		printf("Problem reading frequence.\n");
	else
		printf("At incident angle %f", (float)v);
	theta0 = v*M_PI/180.;

	return;
}
