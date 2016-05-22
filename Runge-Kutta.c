	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "functions.c"

	double *hq(double *q, double h, double tmp, double a, int cnt, double *ans)
	{ 	// реализация веторной операции ans = q + h * tmp / 2 * a, где ans и q - векторы
		int i;
		for (i = 0; i < cnt; i++)
		{
			ans[i] = q[i] + h * tmp / (2 * a);
		}
		return ans;
	}

	double *iter(double *y0, double x0, double a, double h, double *y, int cnt)
	{ // итерация метода Рунге-Кутта второго порядка
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
	{	// реализация веторной операции ans = a + h * k / 2, где ans и a - векторы
		int i;
		for (i = 0; i < cnt; i++)
		{
			ans[i] = a[i] + h * k / 2;
		}
		return ans;
	}

	double *hm(double *a, int cnt, double h, double k, double *ans)
	{
		// реализация веторной операции ans = a + h * k, где ans и a - векторы
		int i;
		for (i = 0; i < cnt; i++)
		{
			ans[i] = a[i] + h * k;
		}
		return ans;
	}

	double *iter2(double *y0, double x0, double h, double *y, int cnt)
	{ // итерация метода Рунге-Кутта четвертого порядка
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

	void Runge_Kutta(double x, double *y0, int cnt, double l, double a, int n)
	{
		int i, j;
		double *y = malloc(cnt * sizeof(double)), h = l / n;
		double *u = malloc(cnt * sizeof(double));
		double *u0 = malloc(cnt * sizeof(double));
		memcpy(u0, y0, cnt * sizeof(double));
		printf("x = %5.2lf, y  = %9lf | y  = %9lf | y  = %9lf\n", x, y0[0], u0[0], proof[0](x));
		//printf("x = %5.2lf, y  = %9lf | y  = %9lf\n", x, y0[0], u0[0]);
		for (i = 1; i < cnt; i++)
			printf("          y%d = %9lf | y%d = %9lf | y%d = %9lf\n",i + 1, y0[i], i + 1, u0[i], i + 1, proof[0](x));
			//printf("          y%d = %9lf | y%d = %9lf\n",i + 1, y0[i], i + 1, u0[i]);
		// y - для метода Рунге-Кутта второго порядка
		// u - для метода Рунге Кутта четвертого порядка 
		for (i = 0; i < n; i++)
		{
			y = iter(y0, x, a, h, y, cnt); //итерация для второго порядка
			u = iter2(u0, x, h, u, cnt); //итерация для четвертого порядка
			x += h;
			printf("x = %5.2lf, ", x);
			printf("y  = %9lf | y  = %9lf | y  =  %9lf\n", y[0], u[0], proof[0](x));
		//	printf("y  = %9lf | y  = %9lf \n", y[0], u[0]);
			for (j = 1; j < cnt; j++)
			{
				printf("           y%d = %9lf | y%d = %9lf | y%d = %9lf\n", j+1, y[j], j+1, u[j], j + 1, proof[0](x));
			//	printf("           y%d = %9lf | y%d = %9lf \n", j+1, y[j], j+1, u[j]);
			}
			memcpy(y0, y, cnt * sizeof(double));
			memcpy(u0, u, cnt * sizeof(double));
		}
		free(y);
		free(u);
		free(u0);
	}

	int main(void)
	{
		double a;
		double x0, *y0, l;
		int cnt, n;
		// считываем все параметры
		scanf("%lf", &x0);
		scanf("%d", &cnt);
		y0 = malloc(cnt * sizeof(double));
		func = malloc(cnt * sizeof(double));
		proof = malloc(cnt * sizeof(double));
		int i;
		for (i = 0; i < cnt; i++)
		{
			scanf("%lf", &y0[i]);
		}
		scanf("%lf", &l);
		scanf("%d", &n);
		scanf("%lf", &a);
		deffunc(); // определяем функции
		//printf("   x         1-го пордяка | 2-го порядка \n");
		printf("   x         1-го пордяка | 2-го порядка   | Аналит.реш.\n");
		Runge_Kutta(x0, y0, cnt, l, a, n);
		free(y0);
		free(func);
		return 0;
	}
