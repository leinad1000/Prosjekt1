#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"
#include "project1.h"

using namespace std;
using namespace arma;

// Opening object that writes to file
ofstream ofile;

void LowerUpper(int n, double *u, double *v, double *x);

// Main program begins
int main(int argc, char *argv[])
{

    // Declaring constants
    int n;
    double average_relative_error,sum,maximum_value;

    char *outfilename;
    // Reading in outputfile. Abort if there are too few command line arguments
    if (argc<=2){
        cout << "Bad usage: " << argv[0] << "read also outputfile and 'n' on the same line"<<endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
        n = atoi(argv[2]);
    }

    // Opening a file for the program
    ofile.open(outfilename);

    // Declaring dynamical arrays
    int *a = new int[n];
    int *b = new int[n];
    int *c = new int[n];
    double *x = new double[n];
    double *f_tilde = new double[n];
    double *v = new double[n];
    double *c_marked = new double[n];
    double *f_marked = new double[n];
    double *u = new double[n];

    double *relative_error = new double[n];
    double *log_ten_of_relative_error = new double[n];

    /* Declaring dynamical matrix
    double **A = new double*[n];
    for(int i = 0; i < n; i++)
        A[i] = new double[n];*/

    // Calling function
    //poissonb(n,a,b,c,x,f_tilde,v,c_marked,f_marked,u);
    //poissonc(n,a,b,c,x,f_tilde,v,c_marked,f_marked,u);
    LowerUpper(n,u,v,x);

    // Writing to file for plotting (open file in MatLab and Import Data and plot in usual way).
    for(int i=0;i<n;i++){
        ofile << setw(15) << setprecision(8) << x[i] << "\t"
              << setw(15) << setprecision(8) << v[i] << "\t"
              << setw(15) << setprecision(8) << u[i] << endl;
    }

    // Question 1d

    // Computing average relative error, for given n
    for(int i=0;i<n;i++){
        log_ten_of_relative_error[i] = log10(fabs((v[i] - u[i])/u[i])); // Computing the 10-logarithm of the relative error for each point
    }

    sum = 0.0;
    for(int i=0;i<n;i++){
        relative_error[i] = (double) pow(10., log_ten_of_relative_error[i]); // Converting back to relative error
        sum += relative_error[i]; // The sum to be used when computing AVERAGE relative error
    }

    average_relative_error = sum/((double) n); // Calculating the average relative error

    // Finding the maximum value of the relative error
    maximum_value = relative_error[0];

    for(int i=1;i<n;i++){
        if(relative_error[i] > relative_error[i-1]){
            maximum_value = relative_error[i];
        }
    }

    // Calculating the 10-logarithm of the step size
    double h,log_h;
    h = 1.0/((double) (n+1));
    log_h = log10(h);

    // Printing out the results
    cout << "log10(h) = " << log_h << endl;
    cout << "Average relative error = " << average_relative_error << endl;
    cout << "Maximum value of relative error = " << maximum_value << endl;
    cout << "n = " << n << endl;

    // End question 1d)

    // Freeing memory
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] x;
    delete [] f_tilde;
    delete [] v;
    delete [] c_marked;
    delete [] f_marked;
    delete [] u;
    delete [] relative_error;
    delete [] log_ten_of_relative_error;

    // Closing output file
    ofile.close();

    // Success
    return 0;
}

// Question e), the LU-decomposition

void LowerUpper(int n,double* u, double*v, double*x){

    // Calculating the step size
    double h;
    h = 1.0/((double) (n+1));
    double hh = h*h;

    // Declaring matrix and vectors
    mat A = zeros<mat>(n,n);
    vec f_tilde(n);
    vec x_points(n);
    vec closed_solution(n);

    // First row of matrix
    A(0,0) = 2;
    A(0,1) = -1;

    // Last row of matrix
    A(n-1,n-2) = -1;
    A(n-1,n-1) = 2;

    // First and last x-value
    x_points(0) = h;
    x_points(n-1) = x_points(0)+(n-1)*h;

    // First and last f-tilde-value
    f_tilde(0) = hh*100*exp(-10*x_points(0));
    f_tilde(n-1) = hh*100*exp(-10*x_points(n-1));

    // Filling up matrix and both vectors completely
    for(int i=1;i<n-1;i++){
            x_points(i) = x_points(i-1) + h;
            f_tilde(i) = hh*100*exp(-10*x_points(i));
            A(i,i-1) = -1;
            A(i,i) = 2;
            A(i,i+1) = -1;
    }

    // Here we start the time-taking
    clock_t start, finish;
    start = clock();

    // Solving by LU-decomposition
    vec solution = solve(A,f_tilde);

    // Here we end the time-taking, FLOPS n*n
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC ); // Calculating time
    cout << "Time LU-decomposition: " << timeused << endl; // Printing out time spent on first part of Thomas algorithm

    // Closed form solution
    double constant;
    constant = 1 - exp(-10);

    for(int i=0;i<n;i++){
        closed_solution(i) = 1 - constant*x_points(i) - exp(-10*x_points(i));
    }

    // Writing the results to vectors that are defined globally
    for(int i=0;i<n;i++){
        v[i] = solution(i);
        u[i] = closed_solution(i);
        x[i] = x_points(i);
    }

}

