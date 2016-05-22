#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "progfunc.c"

void solut(double ct1, double ct2, double ct, double dt1, double dt2, double dt, double start, double end, int n)
{
	double *s, *w, *y;
	s = malloc((n + 1) * sizeof(double));
	w = malloc((n + 1) * sizeof(double));
	y = malloc((n + 1) * sizeof(double));
	double h = (end - start) / n, x = start;
	int i;
	double c, b, a;
	//находим начальные значения прогоночных коэффициенов
	s[1] = -ct2 / (ct1 * h - ct2);
	w[1] = ct * h / (ct1 * h - ct2);
	for (i = 1; i < n; i++)
	{// Прямой ход метода прогонки. Вычисляются прогоночные коэффициенты
		x += h;
		a = (1.0 / (h * h)) - (p(x) / (2.0 * h));
		b = (2.0 / (h * h)) - q(x);
		c = (1.0 / (h * h)) + (p(x) / (2.0 * h));
		s[i + 1] = c / (b - a * s[i]);
		w[i + 1] = (w[i] * a - f(x)) / (b - a * s[i]);
	}
	
	y[n] = (dt2 * w[n] + dt * h) / (dt2 * (1.0 - s[n]) + dt1 * h);
	for (i = n ; i > 0; i--)
	{
		// Обратный ход метода прогонки. По прогоночным коэффициентам, используя
		//рекурентное соотношение находится вектор y
		y[i - 1] = y[i] * s[i] + w[i];
	}
	x = start;
	for (i = 0;  i <= n; i++)
	{
		printf("y(%lf) = %lf\n", x, y[i]);
		x += h; 
	}
}

int main (void)
{
	double ct1, ct2, ct, dt1, dt2, dt, a, b;
	int n;
	// считываем параметры
	scanf("%lf %lf %lf", &ct1, &ct2, &ct);
	scanf("%lf %lf %lf", &dt1, &dt2, &dt);
	scanf("%lf %lf %d", &a, &b, &n);
	deffunc(); // определяем функции
	solut(ct1, ct2, ct, dt1, dt2, dt, a, b, n); // решаем уравненим
	return 0;
}
