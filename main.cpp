#include "diffscatt.h"
#include <iostream>
#include <stdio.h>
#include <vector>

using namespace std;

// Number of points to evaluate
int 		N = 179;
double theta0 = 0;
double angMin = -89;
double angMax = 89;

// Function declaration
void readInput();


int main ( int argc, char *argv[] )
{
	readInput();

	// Iterating
	vector<double> I11;
	vector<double> I22;
	vector<double> I31;
	iterate11(&I11, N, theta0, angMin, angMax);
	iterate22(&I22, N, theta0, angMin, angMax);
	iterate31(&I31, N, theta0, angMin, angMax);

	// Writing output
	FILE * File;
	File = fopen ("out.dat","w");
	for(int i=0; i<N; ++i)
	{
		double theta = angMin+i*(angMax-angMin)/(N-1);
		fprintf(File, "%e\t%e\t%e\t%e\n", theta, I11[i], I22[i], I31[i]);
	}
	fclose (File);

	return 0;
}

void readInput ()
{
	 FILE *input = fopen("par.in", "r");
     if ( input == NULL )
     {
         printf("Unable to read inputfile \'par.in\'");
         return;
     }
 
    char line[80];
    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%d", &N) )
        printf("Problem reading number of points to evaluate.\n");
    else
        printf("Will evaluate %d points.\n", N);
 
    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &theta0) )
        printf("Problem reading incident angle.\n");
    else
		printf("At incident angle %.1f\n", theta0);

    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%d", &corr) )
        printf("Problem reading correlation function.\n");
    else
        cout << "Using correlationfunction " << corr << endl;
 
    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &a) )
        printf("Problem reading transverse correlation length.\n");
    else
        cout << "Using a = " <<  a << endl;

    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &km) )
        printf("Problem reading power spectrum parameter k-.\n");
    else
        cout << "Using k- = " << km << endl;

    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &kp) )
        printf("Problem reading power spectrum parameter k+.\n");
    else
        cout << "Using k+ = " << kp << endl;

	double Re, Im;
    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf+i%lf", &Re, &Im) )
        printf("Problem reading epsilon.\n");
    else
	{
		cout << "Using epsilon = " << Re << " + " << Im << "i" << endl;
		eps = complex<double>(Re,Im);
	}

    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &sigma) )
        printf("Problem reading rms height sigma.\n");
    else
        cout << "Using sigma = " << sigma << endl;

    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &angMin) )
        printf("Problem reading rms height Minimum angle.\n");
    else
        cout << "Using angMin = " << angMin << endl;

    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &angMax) )
        printf("Problem reading rms height Maximum angle.\n");
    else
        cout << "Using angMax = " << angMax << endl;

    return;
}
