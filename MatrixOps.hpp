/*
 * MatrixOps.h
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#ifndef MATRIXOPS_H_
#define MATRIXOPS_H_

#include "User_types.hpp"

void norm_eig_vec(double** EigVectors, int N, d_data domain_data);
void set_matrix(double** A, d_data domain_data, p_params physical_params, s_data *solver_data);
void MatTranspose(double** A, int rows, int columns, double** result);
void VectorNormalization(double *wp, int sizeVec);
void MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2, double** result);
void EigenDecomposition(double** A, int N, double** EigVectors, double* EigValues, int iterations);

#endif /* MATRIXOPS_H_ */
