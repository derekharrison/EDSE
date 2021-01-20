/*
 * main.cpp
 *
 *  Created on: Jan 19, 2021
 *      Author: d-w-h
 *
 *      Implementation of the QR algorithm for computing eigen
 *      vectors and eigen values of a given real symmetric matrix
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "MatrixOps.hpp"
#include "Memory.hpp"

int main(int argc, char* argv[]) {

    double **A, **EigVectors, **Vec, *EigValues, L, h, m, omega;
    int N, iterations;

    /* Parameters */
    N = 50;
    iterations = 5000;
    L = 10.0;
    h = 1.0;
    m = 1.0;
    omega = 1.0;

    /* Some calculations */
    double dx = L/N;
    double alpha = h*h/2/m;
    double *x_p = new double[N];
    double *beta = new double[N];

    /* Compute x_p */
    for(int i = 0; i < N; ++i) {
        x_p[i] = i*dx + 0.5*dx - 0.5*L;
    }

    /* Compute beta */
    for(int i = 0; i < N; ++i) {
        beta[i] = 0.5*m*omega*omega*x_p[i]*x_p[i];
    }

    /* Allocate memory for data */
    A = matrix2D(N, N);
    EigVectors = matrix2D(N, N);
    Vec = matrix2D(N, N);
    EigValues = matrix1D(N);

    /* Initialize and set A */
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j) {
            A[i][j] = 0.0;
        }

    /* First node */
    A[0][0] = 3*alpha/(dx*dx) + beta[0];
    A[0][1] = -alpha/(dx*dx);
    /* Central nodes */
    for(int i = 1; i < N - 1; ++i) {
        A[i][i-1] = -alpha/(dx*dx);
        A[i][i] = 2*alpha/(dx*dx) + beta[i];
        A[i][i+1] = -alpha/(dx*dx);
    }
    /* Last node */
    A[N-1][N-1] = 3*alpha/(dx*dx) + beta[N-1];
    A[N-1][N-2] = -alpha/(dx*dx);

    EigenDecomposition(A, N, EigVectors, EigValues, iterations);

    /* Print results */
    for(int i = 0; i < N; ++i) {
        printf("EigValues[%i]: %f\n", i, EigValues[i]);
    }

    /* Verify computation of eigen vectors */
    MatMult(A, N, N, EigVectors, N, N, Vec);

    printf("Verify computation\n");
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            printf("ratio: %f, ", Vec[i][j]/(EigVectors[i][j] + 1e-20));
        }
        printf("\n");
    }

    /* Normalize eigen vectors */
    double *sum_vec = new double[N];
    for(int j = 0; j < N; ++j) {
        sum_vec[j] = 0.0;
        for(int i = 0; i < N; ++i) {
            sum_vec[j] += EigVectors[i][j]*EigVectors[i][j]*dx;
        }
    }

    for(int j = 0; j < N; ++j) {
        for(int i = 0; i < N; ++i) {
            EigVectors[i][j] = EigVectors[i][j]/sqrt(sum_vec[j]);
        }
    }

    /* Compare with analytical solutions of first three energy levels */
    double *vec_ana = new double[N];
    double *vec_ana_1 = new double[N];
    double *vec_ana_2 = new double[N];
    for(int i = 0; i < N; ++i) {
        vec_ana[i] = pow(m*omega/(M_PI*h),0.25)*exp(-m*omega*x_p[i]*x_p[i]/(2*h));
        vec_ana_1[i] = -pow(m*omega/(M_PI*h),0.25)*sqrt(2*m*omega/h)*x_p[i]*exp(-m*omega*x_p[i]*x_p[i]/(2*h));
        vec_ana_2[i] = -(1-2*m*omega*x_p[i]*x_p[i]/h)*exp(-m*omega*x_p[i]*x_p[i]/(2*h));
    }

    /* Normalize analytical solution level 2 */
    double sum_ana = 0.0;
    for(int i = 0; i < N; ++i) {
        sum_ana += vec_ana_2[i]*vec_ana_2[i]*dx;
    }

    for(int i = 0; i < N; ++i) {
        vec_ana_2[i] = vec_ana_2[i]/sqrt(sum_ana);
    }

    printf("Compare with analytical solution\n");
    for(int i = 0; i < N; ++i) {
        printf("ana_1[%i]: %f, num_1[%i]: %f, ana_2[%i]: %f, num_2[%i]: %f, ana_3[%i]: %f, num_3[%i]: %f\n",
                i, vec_ana[i], i, EigVectors[i][N-1], i, vec_ana_1[i], i, EigVectors[i][N-2], i, vec_ana_2[i], i, EigVectors[i][N-3]);
    }

    /* Export data */
    std::ofstream file;
    std::string file_name = "data.txt";
    file.open(file_name);
    for(int i = 0; i < N; ++i) {
    file << x_p[i] << " "
         << vec_ana[i] << " "
         << EigVectors[i][N-1] << " "
         << vec_ana_1[i] << " "
         << EigVectors[i][N-2] << " "
         << vec_ana_2[i] << " "
         << EigVectors[i][N-3] << "\n";
    }
    file.close();

    /* Free allocated data */
    free2D(A, N);
    free2D(EigVectors, N);
    free2D(Vec, N);
    delete [] EigValues;

    return 0;
}
