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

class DiffScatt
{
  private:
	double power(double k);
	complex<double> alpha (double k);
	complex<double> alpha0 (double k);
	complex<double> G0 (double k);
	complex<double> A1 (double q, double k);
	complex<double> A2 (double q, double p, double k);
	complex<double> A3 (double q, double p, double r, double k);
	double integrandA311Re (double p, void * params);
	double integrandA311Im (double p, void * params);
	complex<double> A311 (double q, double k);
	double integrandT22 (double p, void * params);
	double Ixx (double q, double k, double Txx);

  public:
	void iterate11(double * result, int n, double = -89, double = 89);
	void iterate22(double * result, int n, double = -89, double = 89);
	void iterate31(double * result, int n, double = -89, double = 89);

	// Surface variables
	double a = 100e-9;
	double sigma = 5e-9;
	complex<double> eps (-7.5, 0.24);
}

#endif
