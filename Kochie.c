#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.c"

typedef double(*func_t)(double x, double *y);
func_t *func;

void deffunc(void)
{
	func[0] = f191;
	func[1] = f192;
}

double *hq(double *q, double h, double tmp, double a, int cnt, double *ans)
{
	int i;
	for (i = 0; i < cnt; i++)
	{
		ans[i] = q[i] + h * tmp / (2 * a);
	}
	return ans;
}

double *iter(double *y0, double x0, double a, double h, double *y, int cnt)
{
	int i;
	double tmp;
	double *hlp = malloc (cnt * sizeof(double));
	for (i = 0; i < cnt; i++)
	{
		tmp = func[i](x0, y0);
		y[i] = y0[i] + ((1 - a) * tmp + a * func[i](x0 + h / (2 * a),  hq(y0, h,tmp, a, cnt, hlp))) * h;
	}
	free(hlp);
	return y;
}

double *hn(double *a, int cnt, double h, double k, double *ans)
{
	int i;
	for (i = 0; i < cnt; i++)
	{
		ans[i] = a[i] + h * k / 2;
	}
	return ans;
}

double *hm(double *a, int cnt, double h, double k, double *ans)
{
	int i;
	for (i = 0; i < cnt; i++)
	{
		ans[i] = a[i] + h * k;
	}
	return ans;
}

double *iter2(double *y0, double x0, double h, double *y, int cnt)
{
	int i;
	double k1,k2,k3,k4;
	double *tmp = malloc (cnt * sizeof(double));
	for (i = 0; i < cnt; i++)
	{
		k1 = func[i](x0, y0);
		k2 = func[i](x0 + h/2, hn(y0, cnt, h, k1, tmp));
		k3 = func[i](x0 + h/2, hn(y0, cnt, h, k2, tmp));
		k4 = func[i](x0 + h, hm(y0, cnt, h, k3, tmp));
		y[i] = y0[i] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
	}
	free(tmp);
	return y;
}

void Runge_Kutta4(double x0, double *y0, int cnt, double l, int n)
{
	int i, j;
	printf("x = %lf ", x0);
	for (i = 0; i < cnt; i++)
		printf("   y%d(%lf) = %lf",i, x0, y0[i]);
	printf("\n");
	double *y = malloc(cnt * sizeof(double)), h = l / n;
	for (i = 0; i < n; i++)
	{
		y = iter2(y0, x0, h, y, cnt);
		x0 += h;
		printf("x = %lf ", x0);
		for (j = 0; j < cnt; j++)
			printf("   y%d(%lf) = %lf", j, x0, y[j]);
		printf("\n");
		memcpy(y0, y, cnt * sizeof(double));
	}
	free(y);
}


void Runge_Kutta(double x0, double *y0, int cnt, double l, double a, int n)
{
	int i, j;
	printf("x = %lf ", x0);
	for (i = 0; i < cnt; i++)
		printf("   y%d = %lf",i, y0[i]);
	printf("\n");
	double *y = malloc(cnt * sizeof(double)), h = l / n;
	for (i = 0; i < n; i++)
	{
		y = iter(y0, x0, a, h, y, cnt);
		x0 += h;
		printf("x = %lf ", x0);
		for (j = 0; j < cnt; j++)
			printf("   y%d = %lf",j, y[j]);
		printf("\n");
		memcpy(y0, y, cnt * sizeof(double));
	}
	free(y);
}

int main(void)
{
	double a;
	double x0, *y0, l;
	int cnt, n;
	scanf("%lf", &x0);
	scanf("%d", &cnt);
	y0 = malloc(cnt * sizeof(double));
	double *y0q = malloc(cnt * sizeof(double));
	func = malloc(cnt * sizeof(double));
	int i;
	for (i = 0; i < cnt; i++)
	{
		scanf("%lf", &y0[i]);
	}
	memcpy(y0q, y0, cnt * sizeof(double));
	scanf("%lf", &l);
	scanf("%d", &n);
	scanf("%lf", &a);
	deffunc();
	puts("Метод Рунге-Кутта второго порядка\n");
	Runge_Kutta(x0, y0, cnt, l, a, n);
	puts("--------");
	puts("Метод Рунге-Кутта четвертого порядка\n");
	Runge_Kutta4(x0, y0q, cnt, l, n);
	return 0;
}
