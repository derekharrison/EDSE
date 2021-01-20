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

    double **A, **EigVectors, **Vec, *EigValues, L;
    int N, iterations;
    bool print_verification_results;

    /* Parameters */
    N = domain_data.N;
    L = domain_data.L;

    iterations = 5000;
    print_verification_results = false;

    /* Some calculations */
    double dx = L/N;

    /* Compute x_p */
    for(int i = 0; i < N; ++i) {
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
                printf("ratio: %f, ", (Vec[i][j]/(EigVectors[i][j] + 1e-20))/EigValues[j]);
            }
            printf("\n");
        }
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
    file << solver_data->x_p[i] << " "
         << EigVectors[i][N-1] << " "
         << EigVectors[i][N-2] << " "
         << EigVectors[i][N-3] << "\n";
    }
    file.close();

    /* Free allocated data */
    free2D(A, N);
    free2D(EigVectors, N);
    free2D(Vec, N);
    delete [] EigValues;
}
