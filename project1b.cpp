#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "project1.h"

using namespace std;
extern ofstream ofile;

void poissonb(int n, int* a, int* b, int* c, double* x,double* f_tilde,double* v,double* c_marked,double* f_marked,double* u){

    // Calculating the step size
    double h;
    h = 1.0/((double) (n+1));
    double hh = h*h;

    // Filling up the f-function and its x-values
    x[0] = h; // First term
    for(int i=1;i<n;i++){
        x[i] = x[i-1] + h;
    }
    for(int i=0;i<n;i++){
        f_tilde[i] = hh*100*exp(-10*x[i]);
    }

    // Filling up arrays of the matrix
    for(int i=1;i<n;i++){
        a[i] = -1;
    }
    for(int i=0;i<(n-1);i++){
        c[i] = -1;
    }
    for(int i=0;i<n;i++){
        b[i] = 2;
    }

    // Here we start the time-taking
    clock_t start, finish;
    start = clock();

    // Filling up the constants in the Thomas algorithm; forward substitution
    c_marked[0] = ((double) c[0])/((double) b[0]); // First term
    f_marked[0] = f_tilde[0]/((double) b[0]); // First term

    for(int i=1;i<n;i++){
        c_marked[i] = ((double) c[i])/(b[i] - c_marked[i-1]*a[i]);
        f_marked[i] = (f_tilde[i] - f_marked[i-1]*a[i])/(b[i] - a[i]*c_marked[i-1]);
    } // FLOPS = 8n

    // Solving The Poisson equation by backward substitution
    v[n-1] = f_marked[n-1]; // Last term
    for(int i=(n-2);i>=0;i--){
        v[i] = f_marked[i] - c_marked[i]*v[i+1];
    }

    // Here we end the time-taking, FLOPS = 10n
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC ); // Calculating time
    cout << "Time method b: " << timeused << endl; // Printing out time spent on first part of Thomas algorithm

    // Closed form solution
    double A;
    A = 1 - exp(-10);

    for(int i=0;i<n;i++){
        u[i] = 1 - A*x[i] - exp(-10*x[i]); // FLOPS = 5n
    }

}
