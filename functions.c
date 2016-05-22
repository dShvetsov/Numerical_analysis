#include <math.h>

typedef double(*func_t)(double , double *);
typedef double(*proof_t)(double );

func_t *func;
proof_t *proof;

double f1(double x, double *y)
{
	return 2 * y[0] + y[1];
}

double f2(double x, double *y)
{
	return 3 * y[0] + 4 * y[1];
}
 
double pr1(double x)
{return -exp(x) + exp(5 * x);}
double pr2(double x)
{return exp(x) + 3 * exp(5 * x); } 
 
void deffunc(void)
{
	func[0] = f1;
	func[1] = f2;
	proof[0] = pr1;
	proof[1] = pr2;
}


double fe2(double x, double *y)
{
	return sin(x) - y[0];
}

double fe3(double x, double *y)
{
	return -y[0] - x * x;
}

double fe4(double x, double *y)
{
	return y[0] - y[0] * x;
}

double fe5(double x, double *y)
{
	return (y[0] - y[0] * y[0]) * x;
}

double fe6(double x, double *y)
{
	return (x - x * x) * y[0];
}

double f11(double x, double *y)
{
	return (y[0] - y[1]) / x; 
}

double f12(double x, double *y)
{
	return (y[0] + y[1]) / x;
}


/*double function(double x, double *y)
{
	return (x - x * x) * y[0];
}
*/
