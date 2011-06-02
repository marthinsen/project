/*
 * =====================================================================================
 *
 *       Filename:  diffScatt.h
 *
 *    Description:  Computes the mean diffusive 
 *
 *        Version:  1.0
 *        Created:  30/05/11 11:20:49
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Eirik Marthinsen 
 *        Company:  NTNU
 *
 * =====================================================================================
 */

#ifndef DIFFSCATT_H
#define DIFFSCATT_H

#include <complex>
#include <vector>
using namespace std;

// Call these to calculate the diffuse scattering intensity
void iterate11 (vector<double>* I11, int n, double theta = 0, double min = -89, double max = 89);
void iterate22 (vector<double>* I22, int n, double theta = 0, double min = -89, double max = 89);
void iterate31 (vector<double>* I31, int n, double theta = 0, double min = -89, double max = 89);

// Surface variables
extern double a, sigma, km, kp, km2, kp2;
extern complex<double> eps;
extern int corr;

// Incoming wave
extern double omega;

// Functions
double power(double k);
complex<double> alpha (double k);
complex<double> alpha0 (double k);
complex<double> G0 (double k);
complex<double> A1 (double q, double k);
complex<double> A2 (double q, double p, double k);
complex<double> A3 (double q, double p, double r, double k);
complex<double> A311 (double q, double k);
double integrandA311Re (double p, void * params);
double integrandA311Im (double p, void * params);
double integrandT22 (double p, void * params);
double Ixx (double q, double k, double Txx);

#endif
