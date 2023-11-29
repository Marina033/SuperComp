#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <cstdlib>
#include <cmath>
using namespace std;  

const int N = 10;
const int NN = (N+1)*(N+1);
const int max_it_num = 1000000;

void initializeVariables(double*& a, double*& b, double*& F, double*& w, double*& r, double& x1, double& x2, double& y1, double& y2, double& h1, double& h2, double& metrnorm, int& num_threads, double*& deltas, char* argv[])
{

    int max_it_num = 1000000;
    a = new double [NN]();  // Allocate memory for array 'a'
    b = new double [NN]();  // Allocate memory for array 'b'
    F = new double [NN]();  // Allocate memory for array 'F'
    w = new double [NN]();  // Allocate memory for array 'w'
    r = new double [NN]();  // Allocate memory for array 'r'
    x1 = -4, x2 = 4, y1 = -1, y2 = 4;   //boundaries of rectangle
    h1 = (x2-x1)/N, h2 = (y2-y1)/N;     //grid steps
    metrnorm = h1 * h2; //norm in discr func space

    num_threads = atoi(argv[1]);  // Get the number of threads from command line argument
    omp_set_num_threads(num_threads);  // Set the number of threads for OpenMP
    deltas = new double [num_threads]();  // Allocate memory for array 'deltas'
}

double& get_a(double* a, int i, int j) {
    return a[i * (N + 1) + j];
}

double& get_b(double* b, int i, int j) {
    return b[i * (N + 1) + j];
}

double& get_w(double* w, int i, int j) {
    return w[i * (N + 1) + j];
}

double& get_F(double* F, int i, int j) {
    return F[i * (N + 1) + j];
}

double& get_r(double* r, int i, int j) {
    return r[i * (N + 1) + j];
}

double F1(double x, double y) {
  if (x >= -3 && x <= 3 && y >= 0 && y <= 3) {
    if (y >= (3.0 / 4.0) * x && y >= (-3.0 / 4.0) * x) {
      return 1;
    }
  }
  return 0;
}

/**
 * Calculates the values of 'a' and 'b' using the given parameters.
 *
 * @param a a pointer to the array of 'a' values
 * @param b a pointer to the array of 'b' values
 * @param x1 the starting value of 'x'
 * @param h1 the step size for 'x'
 * @param y1 the starting value of 'y'
 * @param h2 the step size for 'y'
 * @param metrnorm the value of metrnorm
 * @param N the number of iterations
 */
void calculateAB(double* a, double* b, double x1, double h1, double y1, double h2, double metrnorm, int N)
{
    #pragma omp parallel for shared(a, b, x1, h1, y1, h2, metrnorm)
    for (int i = 1; i <= N; i++)
    {
        double x = x1 + (i - 0.5) * h1;
        double yl = (3.0 / 4.0) * x;
        double y = y1 + (i - 0.5) * h2;
        double xl = (1.25 * y);
        for (int j = 1; j <= N; j++)
        {
            double ya = y1 + (j - 0.5) * h2;
            double yb = ya + h2;
            if (yb >= -yl)
                yb = min(yl, yb);
            else
            {
                get_a(a, i, j) = 1.0 / metrnorm;
                get_b(b, i, j) = 1.0 / metrnorm;
                continue;
            }
            if (ya <= yl)
                ya = max(-yl, ya);
            else
            {
                get_a(a, i, j) = 1.0 / metrnorm;
                get_b(b, i, j) = 1.0 / metrnorm;
                continue;
            }
            double l = (yb - ya) / h2;
            get_a(a, i, j) = l + (1 - l) / metrnorm;

            double xa = x1 + (j - 0.5) * h1;
            double xb = xa + h1;
            if (xb >= -xl)
                xb = min(xl, xb);
            else
            {
                get_b(b, i, j) = 1.0 / metrnorm;
                continue;
            }
            if (xa <= xl)
                xa = max(-xl, xa);
            else
            {
                get_b(b, i, j) = 1.0 / metrnorm;
                continue;
            }
            l = (xb - xa) / h1;
            get_b(b, i, j) = l + (1 - l) / metrnorm;
        }
    }
}

/**
 * Saves the results to a file.
 *
 * @param w An array of type double that contains the results to be saved.
 * @param N The size of the array.
 *
 * @return void
 *
 * @throws None
 */
void saveResultsToFile(double* w, int N)
{
    ofstream rslt;
    rslt.open("output.txt");
    for (int j = 0; j <= N; j++)
    {
        rslt << get_w(w, 0, j);
        for (int i = 1; i <= N; i++)
            rslt << ",\t" << get_w(w, i, j);
        rslt << endl;
    }
    rslt.close();
}

/**
 * Initializes variables and performs calculations for a numerical simulation.
 *
 * @param argc the number of command-line arguments
 * @param argv an array of command-line arguments
 *
 * @return an integer representing the exit status of the program
 *
 * @throws None
 */
int main(int argc,char *argv[])
{
    double* a;
    double* b;
    double* F;
    double* w;
    double* r;
    double x1, x2, y1, y2;
    double h1, h2;
    double metrnorm;
    int num_threads;
    double* deltas;

    initializeVariables(a, b, F, w, r, x1, x2, y1, y2, h1, h2, metrnorm, num_threads, deltas, argv);

    
    double time = omp_get_wtime();

    calculateAB(a, b, x1, h1, y1, h2, metrnorm, N);

    #pragma omp parallel for shared(F,x1,h1,y1,h2)
    for (int i=1; i<N; i++)
    {
        double xa=x1+(i-0.5)*h1;
        double xb=xa+h1;
        for (int j=1; j<N; j++)
        {
            double ya=y1+(j-0.5)*h2;
            double yb=ya+h2;
            int n=100;
            double dx=(xb-xa)/n;
            double S=0;
            for (int k=0; k<n; k++)
            {
                double x=xa+(k+0.5)*dx;
                double yl = (3.0 / 4.0) * x;
                S += dx*max(0.0,min(yb,yl)-max(ya,-yl));
            }
            get_F(F,i,j)=S/h1/h2;
        }
    }

    int it=0;
    double dispr=0;
    for (it=0;it<max_it_num; it++)
    {
        #pragma omp parallel for shared(a,b,F,w,r,h1,h2)
        for (int i=1; i<N; i++)
        {
            for (int j=1; j<N; j++)
            {
                get_r(r,i,j) =  get_F(F,i,j);
                get_r(r,i,j) += (get_a(a,i+1,j)*get_w(w,i+1,j) + get_a(a,i,j)*get_w(w,i-1,j))/h1/h1;
                get_r(r,i,j) += (get_b(b,i,j+1)*get_w(w,i,j+1) + get_b(b,i,j)*get_w(w,i,j-1))/h2/h2;
            }
        }


        for (int i=0; i<num_threads; i++) deltas[i]=0;
        #pragma omp parallel for shared(a,b,r,w,h1,h2,deltas)
        for (int i=1; i<N; i++)
        {
            for (int j=1; j<N; j++)
            {
                double wo=get_w(w,i,j);
                get_w(w,i,j) = get_r(r,i,j)/((get_a(a,i+1,j) + get_a(a,i,j))/h1/h1 + (get_b(b,i,j+1) + get_b(b,i,j))/h2/h2);
                double dw=abs(get_w(w,i,j)-wo);
                int id = omp_get_thread_num();
                if ( deltas[id] < dw ) deltas[id] = dw;
            }
        }

        dispr=0;
        for (int i=0; i<num_threads; i++) dispr=max(dispr,deltas[i]); 
        if (dispr<1e-6) break;    
    
}
    cout << "Duration time = " << omp_get_wtime() - time << endl;
    cout << "Iterations = " << it << ", accurancy = " << dispr << endl;
    cout << "Number of threads = " << omp_get_max_threads() << endl;

    ofstream rslt;
    rslt.open ("result.txt");
    for (int j=0; j<=N; j++)
    {
        rslt << w[0,j];
        for (int i=1; i<=N; i++)
            rslt << ",\t" << w[i,j];
        rslt << endl;
    }
    rslt.close();

 
    delete[] a;
    delete[] b;
    delete[] F;
    delete[] w;
    delete[] r;
}