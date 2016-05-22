#include <stdio.h>
#include <stdlib.h>

double **read_matrix(int *n)
{
	int i, j;
	scanf("%d\n", n);
	double **a = calloc (*n, sizeof(double *));
	for (i = 0; i < *n; i++)
	{
		a[i] = calloc (*n + 1, sizeof(double));
		for (j = 0;  j < *n + 1; j++)
		scanf("%lf", &a[i][j]);
	}
	return a;
}

void free_matrix(double **matrix, int n)
{
	int i, j;
	for (i = 0; i < n; i++)
		free(matrix[i]);
	free(matrix);
}

void pr(double **a, int n)
{
	int i, j;
	printf("%d\n", n);
	for (i = 0; i < n; i++)
	{	for (j = 0; j < n + 1; j++)
			printf("%lf ", a[i][j]);
		printf("\n");
	}
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

int main(void)
{
	int n;
	double **a, **b;
	a = read_matrix(&n);
	b = sqr_mtrx(a, n);
	free_matrix(a, n);
	pr(b, n);
	free_matrix(b, n);
	return 0;
}
