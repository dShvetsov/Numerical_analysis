#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#define EPSILON 0.0001
#define qM 1.001 - 2 * M * 0.001

unsigned int count;

enum {
	M = 2,
	N = 40
};

void free_matrix(double **matrix, int n)
{
	int i;
	for (i = 0; i < n; i++)
		free(matrix[i]);
	free(matrix);
}

double **sqr_mtrx(double **a, int n)
{
	int i, j, k;
	double **b = calloc (n, sizeof(double *));
	for (i = 0; i < n; i++)
	{
		b[i] = calloc (n + 1, sizeof(double));
		for (j = 0; j < n; j++)
			b[i][j] = 0;
		b[i][n] = a[i][n];
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				b[i][j] +=a [i][k] * a[k][j];
	return b;
}

double **creat_matrix(int *n, double x)
{
	/*create matrix from of the fomulf*/
	int i, j;
	*n = N;
	double **matrix = calloc (N, sizeof(double *));
	for (i = 0; i < N; i++)
		matrix[i] = calloc (N + 1, sizeof(double));
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			matrix[i][j] = i == j ? pow(qM - 1, i + j) : pow(qM, i + j) + 0.1 * (j + i);
	for (i = 0; i < N; i++)
		matrix[i][N] = fabs(x - (double)N / 10.0) * (i + 1) * sin(x);
	double **hlp = sqr_mtrx(matrix, N);
	free_matrix(matrix, *n);
	return hlp;
}

double **read_matrix(int *n, FILE *f)
{
	int i, j;
	fscanf(f,"%d\n", n);
	double **a = calloc (*n, sizeof(double *));
	for (i = 0; i < *n; i++)
	{
		a[i] = calloc (*n + 1, sizeof(double));
		for (j = 0;  j < *n + 1; j++)
		fscanf(f, "%lf", &a[i][j]);
	}
	return a;
}


double norm(double *a, int n)
{
	/*calculate norm of the vector*/
	int i;
	double ans = 0;
	for (i = 0; i < n; i++)
	{
		ans += a[i] * a[i];
	}
	return sqrt(ans);
}

int accury(double **matrix, double *x, int n)
{
	/* check achieved accury or not*/
	double *dscrp = calloc(n, sizeof(double));
	double tmp = 0;
	int i, j;
	for (i = 0; i < n; i++)
	{
		tmp = 0;
		for (j = 0; j < n; j++)
		{
			tmp += matrix[i][j] * x[j];
		}
		dscrp[i] = tmp - matrix[i][n]; 
	}
	return (norm(dscrp, n) < EPSILON);
}

double *iteration(double **a, double *x, double w, int n, double *x_next)
{
	/* calculate next vector */
	count ++;
	int i, j;
	double tmp1 = 0, tmp2 = 0;
	for (i = 0; i < n; i++)
	{
		tmp1 = 0;
		tmp2 = 0;
		for (j = 0; j <= i - 1; j++)
			tmp1 += a[i][j] * x_next[j];
		for (j = i; j < n; j++)
			tmp2 += a[i][j] * x[j];
		x_next[i] = x[i] + (w / a[i][i]) * (a[i][n]  - tmp1 - tmp2);
	}
	return x_next;
}
void prin(double *x, int n)
{
	int i;
	for (i = 0; i < n; i++)
		printf("%.3lf", x[i]);
	printf("\n");
}

double *upper_relaxation(double **matrix, int n, double w)
{
	/*method of upper_relaxation*/
	int i;
	double *x, *x_next;
	x = calloc(n, sizeof(double));
	x_next = calloc(n, sizeof(double));
	for (i = 0; i < n; i++)
		x[i] = x_next[i] = 0;
	while(!accury(matrix, x, n))
	{
		x_next = iteration(matrix, x, w, n, x_next);
		memcpy(x, x_next, n * sizeof(double));
		//	prin(x, n);
	}
	free(x_next);
	return x;
}

int main(int argc,  char **argv)
{
	FILE *f;
	int n, i, flag;
	double **matrix, *ans;
	double w, x;
	sscanf(argv[2], "%lf", &w);
	sscanf(argv[1], "%d", &flag);
	if (flag == 1)
	{
		f = fopen(argv[3], "r");
		matrix = read_matrix(&n, f);
	}
	else
	{
		sscanf(argv[3],"%lf", &x); 
		matrix = creat_matrix(&n, x);
	}
	puts("ready");
	ans = upper_relaxation(matrix, n, w);
	for (i = 0; i < n; i++)
		printf("X%d = %.3lf; ", i + 1, ans[i]);
	printf("\n");
	free(ans);
	double max_w = 0, min_w = 0;
	unsigned int min = UINT_MAX, max = 0;
	for (w = 0.1; w < 1.95; w += 0.1)
	{
		count = 0;
		free(upper_relaxation(matrix, n, w));
		if (count > max)
		{
			max = count;
			max_w = w;
		}
		if (count < min)
			{
				min = count;
				min_w = w;
			}
		printf("Done %lf\n", w);
	}
	printf("скорость сходимости является максимальной при w = %lf (%d итераций)\n", min_w, min);
	printf("скорость сходимости является минимальной при w = %lf (%d итераций)\n", max_w, max);
	free_matrix(matrix, n);
	return 0; 
}
