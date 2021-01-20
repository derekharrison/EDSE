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
#include <time.h>
#include <math.h>
#include "MatrixOps.hpp"
#include "Memory.hpp"
#include "Solver.hpp"
#include "User_types.hpp"

double V(p_params physical_params, double x) {
    /* Potential */
    double f = 0.5*physical_params.m*physical_params.omega*
                   physical_params.omega*x*x;

    return f;
}

int main(int argc, char* argv[]) {

    d_data domain_data = {0};
    p_params physical_params = {0};
    s_data solver_data = {0};
    clock_t t_start = clock();

    /* Parameters */
    domain_data.N = 50;          //Number of nodes
    domain_data.L = 10.0;        //Length of domain

    physical_params.h = 1.0;     //Constant
    physical_params.m = 1.0;     //Particle mass
    physical_params.omega = 1.0; //Frequency

    /* Allocate memory for solver results */
    solver_data.EigVectors = matrix2D(domain_data.N, domain_data.N);
    solver_data.EigValues = matrix1D(domain_data.N);
    solver_data.x_p = matrix1D(domain_data.N);

    /* Execute solver */
    solver(domain_data, physical_params, &solver_data);

    /* Print results */
    printf("Eigen values:\n");
    for(int i = 0; i < domain_data.N; ++i) {
        printf("EigValues[%i]: %f\n", i, solver_data.EigValues[i]);
    }

    printf("Eigen vectors:\n");
    for(int i = 0; i < domain_data.N; ++i) {
        for(int j = 0; j < domain_data.N; ++j) {
            printf("[%i][%i]: %f, ", i, j, solver_data.EigVectors[i][j]);
        }
        printf("\n");
    }

    clock_t t_end = clock();
    double time_execution = double (t_end - t_start) / CLOCKS_PER_SEC;
    printf("execution time: %f\n", time_execution);

    return 0;
}
