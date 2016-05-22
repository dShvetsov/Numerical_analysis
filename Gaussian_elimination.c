#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPSILON 0.0000001
#define qM 1.001 - 2 * M * 0.001

typedef double **matrix;
typedef double * vector;


enum {
	M = 2,
	N = 40
};
	
enum { 
	EAS = 0, //Gaussian elimination without pivoting
	PVT = 1  //Gaussuan elimination with pivoting
};

typedef int (*find_func)(matrix, int, int);
typedef void(*stdtion_func)(matrix, int, int);
typedef void(*subtract_func)(matrix, int, int, int);
typedef void(*calculate_func)(matrix, int, int *, vector, int);

matrix read_matrix(int *n, FILE *f)
{
	int i, j;
	fscanf(f,"%d\n", n);
	matrix a = calloc (*n, sizeof(vector));
	for (i = 0; i < *n; i++)
	{
		a[i] = calloc (*n + 1, sizeof(double));
		for (j = 0;  j < *n + 1; j++)
			fscanf(f,"%lf", &a[i][j]);
	}
	return a;
}

matrix conversion_matrix(matrix src, int n)
{
	/*conversion matrix from a[string][column] to a[column][string]*/
	int i, j;
	matrix mtrx = calloc(n + 1, sizeof(vector));
	for (i = 0; i < n + 1; i++)
		mtrx[i] = calloc(n, sizeof(double));
	for (i = 0; i < n; i++)
		for (j = 0; j < n + 1; j++)
			mtrx[j][i] = src[i][j];
	return mtrx;
}

matrix creat_matrix(int *n, double x)
{
	/*create matrix from of the formula*/
	int i, j;
	*n = N;
	matrix mtrx = calloc (N, sizeof(vector));
	for (i = 0; i < N; i++)
		mtrx[i] = calloc (N + 1, sizeof(double));
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			mtrx[i][j] = i == j ? pow(qM - 1, i + j) : pow(qM, i + j) + 0.1 * (j - i);
	for (i = 0; i < N; i++)
		mtrx[i][N] = fabs(x - (double)N / 10.0) * (i + 1) * sin(x);
	return mtrx;
}

matrix copy_matrix(matrix src, int n, int m)
{
	int i;
	matrix dst = calloc(n, sizeof(vector));
	for (i = 0; i < n; i++)
	{
		dst[i] = calloc(m, sizeof(double));
		memcpy(dst[i], src[i], m * sizeof(double));
	}
	return dst;
}

matrix attached_matrix(matrix mtrx, int n)
{
	/*creat attached matrix : [A|I]*/
	int i, j;
	matrix mtrx2 = calloc(n, sizeof(vector));
	for (i = 0; i < n; i++)
	{
		mtrx2[i] = calloc(2 * n, sizeof(double));
		memcpy(mtrx2[i], mtrx[i], n * sizeof(double)); 
		for (j = n; j < 2 * n; j++)
			mtrx2[i][j] = (double)(i == j - n);
	}
	return mtrx2;
}

void free_matrix(matrix mtrx, int n)
{
	int i;
	for (i = 0; i < n; i++)
		free(mtrx[i]);
	free(mtrx);
}

int find_max(matrix mtrx, int n, int numstr)
{
	/*find maximal element in numstr'st string */
	/*need for Gaussian elimination with pivoting */
	/*matrix[column][string] */
	int i, tmp = numstr;
	double max = fabs(mtrx[numstr][numstr]);
	for (i = numstr; i < n; i++)
	{
		if (max < fabs(mtrx[i][numstr]))
		{
			max = fabs(mtrx[i][numstr]);
			tmp = i;
		}
	}
	if (max < EPSILON)
		puts("matrix is degenerate");
	return tmp;
}

int find_not_null_clm(matrix mtrx, int n, int num)
{
	/*find not null element in num'st column*/
	int i;
	for (i = num; i < n; i++)
	{
		if (fabs(mtrx[num][i]) > EPSILON)
			return i;
	}
	puts("matrix is degenerate");
	return -1;
}

int my_swap(matrix mtrx, int *form, int new, int old)
{
	/*swap new'st and old'st columns (or strings) in matrix*/
	/*if necessary describe swap variables*/
	/*returned value: 1 if wasn't swap, -1 if was swap*/
	if (new == old)
		return 1;
	vector hlp = mtrx[new];
	mtrx[new] = mtrx[old];
	mtrx[old] = hlp;
	if( form != NULL)
	{
		int hlp2 = form[new];
		form[new] = form[old];
		form[old] = hlp2;
	}
	return -1;
}

void stdtion_of_string(matrix mtrx, int numstr, int n)
{
	/*divide by the dioganal element of string*/
	/*for Gaussian elimination with pivoting*/
	/*matrix[column][string]*/
	double a = mtrx[numstr][numstr];
	int i;
	mtrx[numstr][numstr] = 1.0;
	for (i = numstr + 1; i < n; i++)
	{
		mtrx[i][numstr] /= a;
	}
}

void stdtion_of_string2(matrix mtrx, int numstr, int n)
{
	/*divide by the dioganal element of string*/
	int i;
	double a = mtrx[numstr][numstr];
	for (i = numstr; i < n; i++)
	{
		mtrx[numstr][i] /= a;
	}
}

void subtract_string(matrix mtrx, int num1, int num2, int n)
{
	/*subtract num1'st string of num2'st string*/
	/*for Gaussian elimination with pivoting*/
	/*matrix[column][string]*/
	double factor = mtrx[num1][num2];
	int i;
	mtrx[num1][num2] = 0.0;
	for (i = num1 + 1; i < n; i++)
	{
		mtrx[i][num2] -= mtrx[i][num1] * factor;
	}
}

void subtract_string2(matrix mtrx, int num1, int num2, int n)
{
	/*subtract num1'st string of num2'st string*/
	double factor = mtrx[num2][num1];
	int i;
	mtrx[num2][num1] = 0.0;
	for (i = num1 + 1; i < n; i++)
	{
		mtrx[num2][i] -= mtrx[num1][i] * factor;
	}
}

void calculate_ans2(matrix mtrx, int numstr, int *a, vector ans, int n)
{
	/*calculate answer for Gaussian elimination*/
	double tmp = mtrx[numstr][n];
	int i;
	for (i = numstr + 1; i < n; i++)
	{
		tmp -= (mtrx[numstr][i] * ans[i]);
	}
	ans[numstr] = tmp;
}

void calculate_ans(matrix mtrx, int numstr, int *form, vector ans, int n)
{
	/*calculate answer for Gaussian elimination with pivoting*/
	/*matrix[column][sting]*/
	double tmp = mtrx[n][numstr];
	int i;
	for (i = numstr + 1; i < n; i++)
	{
		tmp -= (mtrx[i][numstr] * ans[form[i]]);
	}
	ans[form[numstr]] = tmp;
}

vector Gauss(matrix src, int n, double *determ, int mode) 
{ 
	matrix mtrx;
	int *form, sign = 1;
	find_func find;
	stdtion_func stdtion;
	subtract_func subtract;
	calculate_func calculate;
	//two methods differ in set of functions
	//choose set of function
	if (mode == PVT)
	{
	/*Gaussian elimination with pivoting*/
	/*matrix is represented as matrix[column][string]*/
		int i;
		form = calloc (n, sizeof(int));
		for (i = 0; i < n; i++)
			form[i] = i;
		find = find_max;
		stdtion = stdtion_of_string;
		subtract = subtract_string;
		calculate = calculate_ans;
		mtrx = conversion_matrix(src, n);
	}
	else
	{
	/*Gaussian elimination*/
		form = NULL;
		find = find_not_null_clm;
		stdtion = stdtion_of_string2;
		subtract = subtract_string2;
		calculate = calculate_ans2;
		mtrx = copy_matrix(src, n, n + 1);
	}
	vector ans;
	double det = 1;
	int i, hlp, j;
	for (i = 0; i < n; i++) //forward stroke Gaussian elimination
	{
		hlp = find(mtrx, n, i); 
		sign *= my_swap(mtrx, form, hlp, i);
		if (determ != NULL)
			det *= mtrx[i][i]; //if necessary calculate deterinate
		stdtion(mtrx, i, n + 1);
		for (j = i + 1; j < n; j++)
		{
			subtract(mtrx, i, j, n + 1);
		}
	}
	ans = calloc(n, sizeof(double));
	for (i = n - 1; i >= 0; i--)
	{
		calculate(mtrx, i, form, ans, n); //reversal Gaussian elimination
	}
	if (determ != NULL)
		*determ = sign * det;
	if (form != NULL)
		free(form);
	mode == PVT ? free_matrix(mtrx, n + 1) : free_matrix(mtrx, n);
	return ans;
}


matrix inverse(matrix src, int n)
{
	/*calculate inverse matrix*/
	matrix mtrx = attached_matrix(src, n);
	int tmp, i, j;
	for (i = 0; i < n; i++)
	{
		tmp = find_not_null_clm(mtrx, n, i);
		my_swap(mtrx, NULL, tmp, i);
		stdtion_of_string2(mtrx, i, 2 * n);
		for (j = i + 1; j < n; j++)
		{
			subtract_string2(mtrx, i, j, 2 * n);
		}
	}
	for (i = n - 2; i >= 0; i--)
	{
		for (j = i + 1; j < n; j++)
		{
			subtract_string2(mtrx, j, i, 2 * n);
		}
	}
	matrix answer = calloc(n, sizeof(vector));
	for (i = 0; i < n; i++)
	{
		answer[i] = calloc(n, sizeof(double));
		for (j = 0; j < n; j++)
		{
			answer[i][j]  = mtrx[i][n + j];
		}
	}
	free_matrix(mtrx, n);
	return answer;
}

void report_answer(vector answ1, vector answ2, double det, matrix mtrx, int n)
{
	int i, j;
	puts("Решение СЛАУ методом Гаусса без выбора главного элемента:");
	for (i = 0; i < n; i++)
		printf("X%d = %.3lf; ", i + 1, answ1[i]); 
	puts("\n");
	puts("Решение СЛАУ методом Гаусса с выбором главного элемента:");
	for (i = 0; i < n; i++)
		printf("X%d = %.3lf; ", i+1, answ2[i]);
	puts("\n");
	printf("Определитель матрицы:\ndet(A) = %.3lf\n", det);
	printf("\n");
	puts("Обратная матрица:");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("%.1lf ", mtrx[i][j]);
		printf("\n");
	}
}

int main (int argc, char **argv)
{
	int n;
	if(argc < 2)
	{
		puts("Too few parameters");
		puts("write type of submission matrix (1 -- from file, 2 -- formula)");
		return 0;
	}
	int submission;
	sscanf(argv[1],"%d", &submission);
	double determ;
	vector answer1, answer2;
	//matrix mtrx, **mtrx2, **reverse, **mtrx3;
	matrix mtrx, reverse;
	double x = 1;
	FILE *f;
	if (submission == 1)
	{
		if (argc >=2)
			{if ((f = fopen(argv[2], "r")) == NULL)
			{
				puts("error with file");
				return -1;
			}}
		else
			if( (f = fopen("input.txt", "r")) == NULL)
			{
				puts("please try again and se file");
				return -1;
			}
		mtrx = read_matrix(&n, f);
	}
	else
	{
		if (argc >= 2)
			sscanf(argv[2], "%lf", &x);
		mtrx = creat_matrix(&n, x);
	}
		answer2 = Gauss(mtrx, n, NULL, EAS);
	answer1 = Gauss(mtrx, n, &determ, PVT);
	reverse = inverse(mtrx, n);
	report_answer(answer1, answer2, determ, reverse, n);
	free_matrix(reverse, n); 
	free_matrix(mtrx, n);
	return 0;
}
