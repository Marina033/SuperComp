#ifndef NUMSOL
#define NUMSOL
#define WIDTH_FIELD 6			
#define PRECISION_AFTER_POINT 1 

#include <iomanip>
#include <iostream>
#include <cmath>
#include <omp.h>

using namespace std;


double Uxy(double x, double y)
{
	return 1 - (x - 1) * (x - 1) - (y - 0.5) * (y - 0.5);
}
double f(double x, double y) 
{
    if (x * x + 4 * y * y < 1) {
        return 1;
    } 
	else {
        return 0;
    }
}

double M1(double y)
{
	return 0;
}
double M2(double y)
{
	return 0;
}
double M3(double x)
{
	return 0;
}
double M4(double x)
{
	return 0;
}

double** MemoryAllocator(int n, int m)
{
	double** Matrix = NULL;

	Matrix = new double* [n];
	for (int i = 0; i < n; i++)
		Matrix[i] = new double[m];

	return Matrix;
}
void MemoryCleaner(double** arr, int n)
{
	for (int i = 0; i < n; i++)
		delete[] arr[i];

	delete[] arr;
	arr = NULL;
}
void ShowSolution(double** V, int n, int m)
{
	for (int j = m; j >= 0; j--)
	{
		for (int i = 0; i <= m; i++)
			printf("%5.3lf   ", V[i][j]);
		printf("\n");
	}
}
void FillRightSide_omp(double** F, int n, int m, double a, double c, double h, double k)
{
    for (int j = 1; j < m; j++)
	#pragma omp parallel for
        for (int i = 1; i < n; i++)
        {
            double Xi, Yj, sum = 0;
            Xi = a + i * h;
            Yj = c + j * k;

            if (j == 1)
                sum += (1 / (k * k)) * M3(Xi);
            else if (j == m - 1)
                sum += (1 / (k * k)) * M4(Xi);
            if (i == 1)
                sum += (1 / (h * h)) * M1(Yj);
            else if (i == n - 1)
                sum += (1 / (h * h)) * M2(Yj);

            F[i][j] = -f(Xi, Yj) - sum;
        }
}
void FillRightSide(double** F, int n, int m, double a, double c, double h, double k)
{
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			double Xi, Yj, sum = 0;
			Xi = a + i * h;
			Yj = c + j * k;

			if (j == 1)
				sum += (1 / (k * k)) * M3(Xi);
			else
				if (j == m - 1)
					sum += (1 / (k * k)) * M4(Xi);
			if (i == 1)
				sum += (1 / (h * h)) * M1(Yj);
			else
				if (i == n - 1)
					sum += (1 / (h * h)) * M2(Yj);

			F[i][j] = -f(Xi, Yj) - sum;
		}
}
/**
 * Fills the given 2D array with values based on the provided parameters.
 *
 * @param V The 2D array to be filled.
 * @param n The number of columns in the array.
 * @param m The number of rows in the array.
 * @param a The starting x-coordinate of the array.
 * @param b The ending x-coordinate of the array.
 * @param c The starting y-coordinate of the array.
 * @param d The ending y-coordinate of the array.
 *
 * @throws None.
 */
void FillStartSolution_omp(double** V, int n, int m, double a, double b, double c, double d)
{
    double h, k;
    h = (b - a) / n;
    k = (d - c) / m;

    #pragma omp parallel for
    for (int j = 0; j <= m; j++)
        for (int i = 0; i <= n; i++)
        {
            if (i == 0 || j == 0 || i == n || j == m)
            {
                double Xi, Yj, sum = 0;
                Xi = a + i * h;
                Yj = c + j * k;
                if (j == 0)
                    V[i][j] = M3(Xi);
                else if (j == m)
                    V[i][j] = M4(Xi);
                if (i == 0)
                    V[i][j] = M1(Yj);
                else if (i == n)
                    V[i][j] = M2(Yj);
            }
            else
                V[i][j] = 0;
        }
}
/**
 * Fills the start solution matrix with values based on the given parameters.
 *
 * @param V The matrix to be filled with the start solution values.
 * @param n The number of columns in the matrix.
 * @param m The number of rows in the matrix.
 * @param a The start value for the x-axis.
 * @param b The end value for the x-axis.
 * @param c The start value for the y-axis.
 * @param d The end value for the y-axis.
 *
 * @throws None.
 */
void FillStartSolution(double** V, int n, int m, double a, double b, double c, double d)
{
	double h, k;
	h = (b - a) / n;
	k = (d - c) / m;

	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			if (i == 0 || j == 0 || i == n || j == m)
			{
				double Xi, Yj, sum = 0;
				Xi = a + i * h;
				Yj = c + j * k;
				if (j == 0)
					V[i][j] = M3(Xi);
				else
					if (j == m)
						V[i][j] = M4(Xi);
				if (i == 0)
					V[i][j] = M1(Yj);
				else
					if (i == n)
						V[i][j] = M2(Yj);
			}
			else
				V[i][j] = 0;
		}
}
void ZeidelsMethod_omp(double** V, int n, int m, double a, double b, double c, double d, double eps, int Nmax, double& epsMax, int& S)
{
    double epsCur = 0;
    double a2, k2, h2;
    double v_old;
    double v_new;

    h2 = -((n / (b - a)) * (n / (b - a)));
    k2 = -((m / (d - c)) * (m / (d - c)));
    a2 = -2 * (h2 + k2);

    while (true)
    {
        epsMax = 0;

        for (int j = 1; j < m; j++)
        {
			#pragma omp parallel for private(v_old, v_new, epsCur) shared(V, epsMax)
            for (int i = 1; i < n; i++)
            {
                double Xi, Yj;
                Xi = a + i * ((b - a) / n);
                Yj = c + j * ((d - c) / m);

                v_old = V[i][j];
                v_new = -(h2 * (V[i + 1][j] + V[i - 1][j]) + k2 * (V[i][j + 1] + V[i][j - 1]));
                v_new = v_new + f(Xi, Yj);
                v_new = v_new / a2;

                epsCur = abs(v_old - v_new);
#pragma omp critical
                {
                    if (epsCur > epsMax)
                        epsMax = epsCur;
                    V[i][j] = v_new;
                }
            }
        }

        if (S == 0 || S == 1)
        {
            cout << " " << S + 1 << " S " << endl;
            ShowSolution(V, n, m);
            cout << endl;
        }
        ++S;

        if ((epsMax < eps) || (S >= Nmax))
            break;
    }
}
/**
 * ZeidelsMethod performs Zeidel's method to solve a system of partial differential equations.
 *
 * @param V The 2D array representing the grid of values.
 * @param n The number of grid points in the x-direction.
 * @param m The number of grid points in the y-direction.
 * @param a The lower bound of the x-domain.
 * @param b The upper bound of the x-domain.
 * @param c The lower bound of the y-domain.
 * @param d The upper bound of the y-domain.
 * @param eps The desired error tolerance.
 * @param Nmax The maximum number of iterations.
 * @param epsMax The maximum error encountered during the iterations.
 * @param S The number of iterations performed.
 *
 * @throws None
 */
void ZeidelsMethod(double** V, int n, int m, double a, double b, double c, double d, double eps, int Nmax, double& epsMax, int& S)
{
	double	epsCur = 0;			
	double	a2, k2, h2;			
	double	v_old;				
	double	v_new;				

	h2 = -((n / (b - a)) * (n / (b - a)));
	k2 = -((m / (d - c)) * (m / (d - c)));
	a2 = -2 * (h2 + k2);

	while (true)
	{
		epsMax = 0;
		for (int j = 1; j < m; j++)
			for (int i = 1; i < n; i++)
			{
				double Xi, Yj;
				Xi = a + i * ((b - a) / n);
				Yj = c + j * ((d - c) / m);

				v_old = V[i][j];
				v_new = -(h2 * (V[i + 1][j] + V[i - 1][j]) + k2 * (V[i][j + 1] + V[i][j - 1]));
				v_new = v_new + f(Xi, Yj);
				v_new = v_new / a2;

				epsCur = abs(v_old - v_new);
				if (epsCur > epsMax)
					epsMax = epsCur;

				V[i][j] = v_new;
			}
		if (S == 0 || S == 1)
		{
			cout << " " << S + 1 << " S " << endl;
			ShowSolution(V, n, m);
			cout << endl;
		}
		++S;

		if ((epsMax < eps) || (S >= Nmax))
			break;
	}
}


double DiscrepancyOfSolution(double** V, int n, int m, double a, double b, double c, double d)
/**
 * Calculates the discrepancy of a solution.
 *
 * @param V A 2D array representing the solution matrix
 * @param n The number of rows in the solution matrix
 * @param m The number of columns in the solution matrix
 * @param a The lower bound of the x-axis range
 * @param b The upper bound of the x-axis range
 * @param c The lower bound of the y-axis range
 * @param d The upper bound of the y-axis range
 *
 * @return The discrepancy of the solution
 *
 * @throws None
 */
{
	double	a2, k2, h2;			
	double  h, k;				
	double** F;					
	double rs = 0;				

	h = (b - a) / n;
	k = (d - c) / m;

	h2 = ((n / (b - a)) * (n / (b - a)));
	k2 = ((m / (d - c)) * (m / (d - c)));
	a2 = -2 * (h2 + k2);

	
	F = MemoryAllocator(n + 1, m + 1);
	FillRightSide(F, n, m, a, c, h, k);

	for (int j = 1; j < m; j++)
	{
		for (int i = 1; i < n; i++)
		{
			double r;
			double mult;

			if (j != 1 && j != m - 1)
			{
				
				if (i != 1 && i != n - 1)
					mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
				else
					if (i == 1)
						mult = k2 * V[i][j - 1] + a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
					else
						if (i == n - 1)
							mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j] + k2 * V[i][j + 1];
			}
			else
				if (j == 1)
				{
					if (i == 1)
						mult = a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
					else
						if (i != n - 1)
							mult = h2 * V[i - 1][j] + a2 * V[i][j] + h2 * V[i + 1][j] + k2 * V[i][j + 1];
						else
							if (i == n - 1)
								mult = h2 * V[i - 1][j] + a2 * V[i][j] + k2 * V[i][j + 1];
				}
				else
					if (j == m - 1)
					{
						if (i == 1)
							mult = k2 * V[i][j - 1] + a2 * V[i][j] + h2 * V[i + 1][j];
						else
							if (i != n - 1)
								mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j] + h2 * V[i + 1][j];
							else
								if (i == n - 1)
									mult = k2 * V[i][j - 1] + h2 * V[i - 1][j] + a2 * V[i][j];
					}

			r = abs(mult - F[i][j]);
			rs += r * r;
		}
	}
	MemoryCleaner(F, n);

	return sqrt(rs);
}

double CheckSolution(double** V, int n, int m, double a, double b, double c, double d)
{
	
	double** U = MemoryAllocator(n + 1, m + 1);
	double h, k;
	double zs = 0;

	h = (b - a) / n;
	k = (d - c) / m;

	
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			double Xi, Yj;
			Xi = a + i * h;
			Yj = c + j * k;

			U[i][j] = Uxy(Xi, Yj);
		}
	
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			double z = abs(U[i][j] - V[i][j]);

			if (z > zs)
				zs = z;
		}
	MemoryCleaner(U, n);

	return zs;
}

#endif