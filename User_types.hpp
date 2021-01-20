/*
 * User_types.hpp
 *
 *  Created on: Jan 20, 2021
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

typedef struct domain_data {
    int N;
    double L;
} d_data;

typedef struct physical_params {
    double h;
    double m;
    double omega;
} p_params;

typedef struct solver_data {
    double** EigVectors;
    double* EigValues;
    double* x_p;
} s_data;

#endif /* USER_TYPES_HPP_ */
