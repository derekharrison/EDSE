/*
 * Solver.cpp
 *
 *  Created on: Jan 20, 2021
 *      Author: d-w-h
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "main.hpp"
#include "MatrixOps.hpp"
#include "Memory.hpp"
#include "User_types.hpp"

void solver(d_data domain_data, p_params physical_params, s_data* solver_data) {

    double **A, **EigVectors, **Vec, *EigValues, L, h, m, omega;
    int N, iterations;
    bool print_verification_results;

    /* Parameters */
    N = domain_data.N;
    L = domain_data.L;

    h = physical_params.h;
    m = physical_params.m;
    omega = physical_params.omega;

    iterations = 5000;
    print_verification_results = false;

    /* Some calculations */
    double dx = L/N;
    double *x_p = new double[N];

    /* Compute x_p */
    for(int i = 0; i < N; ++i) {
        x_p[i] = i*dx + 0.5*dx - 0.5*L;
        solver_data->x_p[i] = i*dx + 0.5*dx - 0.5*L;
    }

    /* Allocate memory for data */
    A = matrix2D(N, N);
    EigVectors = matrix2D(N, N);
    Vec = matrix2D(N, N);
    EigValues = matrix1D(N);

    /* Initialize and set A */
    set_matrix(A, domain_data, physical_params, solver_data);

    /* Compute eigen values and vectors */
    EigenDecomposition(A, N, EigVectors, EigValues, iterations);

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

    if(print_verification_results) {
        /* Verify computation of eigen vectors */
        MatMult(A, N, N, EigVectors, N, N, Vec);

        printf("Verify computation\n");
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                printf("ratio: %f, ", Vec[i][j]/(EigVectors[i][j] + 1e-20));
            }
            printf("\n");
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

    /* Set results */
    for(int i = 0; i < N; ++i) {
        solver_data->EigValues[i] = EigValues[i];
        for(int j = 0; j < N; ++j) {
            solver_data->EigVectors[i][j] = EigVectors[i][j];
        }
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
}
