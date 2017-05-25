#define _MSL_H

#ifndef __IOSTREAM__
#include <iostream>
#endif

#ifndef __CASSERT__
#include <cassert>
#endif

#ifndef __CSTDIO__
#include <cstdio>
#endif


#ifndef __CMATH__
#include <cmath>
#endif

#ifndef __FSTREAM__
#include <fstream>
#endif

#ifndef __STRING__
#include <string>
#endif

#ifndef __SGI_STL_VECTOR
#include <vector>
#endif

#ifndef __SGI_STL_ALGORITHM
#include <algorithm>
#endif 

#define Pi 3.14159265
#define pi 3.14156

int int_max(int x, int y);
int int_min(int x, int y);
double double_max(double x, double y);
double double_min(double x, double y);
double square(double x);
double sqr(double x);
double Z(double x);
double P(double x);
double Q(double x);
double I(double x, double a, double b);
double A(double t, int nu);
double t_prob(double t, int nu);
double Gamma(int n);
double logGamma(int n);
double log_factorial(int n);
double left_bin_p_value(int x, double p, int n);
double poiss_p_value(int x, double lambda);
double right_bin_p_value(int x, double p, int n);
double sqr(double x);
double logBesselK_semi_int(double nu, double z);
double logBesselK(int n, double z);
double log_poiss_pr(int x, double lambda);
