/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Compute the mean diffusive intesity of scattered light on rough 
 *    				surface.
 *
 *        Version:  1.0
 *        Created:  31/05/11 12:10:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eirik Marthinsen (eirikma@stud.ntnu.no), 
 *        Company:  NTNU, Trondheim
 *
 * =====================================================================================
 */

#include "diffscatt.h"
#include <iostream>
#include <stdio.h>
#include <vector>

using namespace std;

// Number of points to evaluate
int 		N = 179;
double theta0 = 40;
double angMin = -89;
double angMax = 89;

// Function declaration
void readInput();

/* 
 */
int main ( int argc, char *argv[] )
{
  	// Read inputfile if any
	readInput();

	// Iterating
	vector<double> I11;
	vector<double> I22;
	vector<double> I31;
	iterate11(&I11, N, theta0);
	iterate22(&I22, N, theta0);
	iterate31(&I31, N, theta0);

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

/*
 * Reads inputfile 'par.in'
 */
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
 
    double v;
    if( !fgets(line, sizeof(line), input ) || !sscanf( line, "%lf", &v) )
        printf("Problem reading incident angle.\n");
    else
	{
		printf("At incident angle %f\n", (float)v);
		theta0 = v*M_PI/180.;
	}

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
    return;
}
