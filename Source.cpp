#include "NumberSolution.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

using namespace std;


//U"xx + U"yy = -f(x, y)
//a <= x <= b, c <= y <= d
//U(a, y) = M1(y), U(b, y) = M2(y)
//U(x, c) = M3(x), U(x, d) = M3(x)

int main(void)
{

	system("________________________");

	int		Nmax = 10000;		
	int		S = 0;				
	double	eps = 0.0000001;	
	double	epsMax = 0;			
	int		n = 20, m = 20;		
	double** V = NULL;			
	double** F = NULL;			
	double	a, b, c, d;			
	int Exit = 1, Show = 0;

	a = -1;
	b = 1;
	c = -0.5;
	d = 0.5;

	Nmax = 0;

	system("cls");
	cout << "starting" << endl;

	while (Nmax <= 0)
	{
		cout << endl << "input max iterations: ";
		cin >> Nmax;
	}
	system("cls");

	V = MemoryAllocator(n + 1, m + 1);
	FillStartSolution_omp(V, n, m, a, b, c, d);
	ZeidelsMethod_omp(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
	cout << "solving" << endl;
	ShowSolution(V, n, m);

	//˜˜˜˜˜˜˜
	cout << endl << "---------------------------------------------" << endl;
	cout << "n,m: (" << n << ", " << m << ")" << endl;
	cout << "bounds for X: [" << a << "; " << b << "]" << endl;
	cout << "bounds for Y: [" << c << "; " << d << "]" << endl;
	cout << "grid step Ox: h = " << (b - a) / n << endl;
	cout << "grid step Oy: k = " << (d - c) / m << endl;
	cout << "epsilon: " << eps << endl;
	cout << "eps_max: " << epsMax << endl;
	cout << "S: " << S << endl;
	cout << " e-solution: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;
	cout << "---------------------------------------------" << endl << endl;

	MemoryCleaner(V, n);
	
	cout << endl;
	return 0;
}