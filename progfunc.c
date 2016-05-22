#include <math.h>

typedef double (*func_t) (double x);
func_t p, q, f, proof;

double p1(double x)
{return 2 * x * x;}

double q1(double x)
{return 1;}

double f1(double x)
{return x;}

void deffunc(void)
{
	p = p1;
	q = q1; 
	f = f1;
}


double pz(double x)
{return 2;}

double qz(double x)
{return 2;}

double fz(double x)
{return x* exp(-x);}

double pr(double x)
{return exp(-x) * (x - sin(x));}



double p2(double x)
{return 0;}

double q2(double x)
{return -1;}

double f2(double x)
{return 2 * x;}

double pr2(double x)
{return  (sinh(x)/sinh(1)) - 2 * x;}

