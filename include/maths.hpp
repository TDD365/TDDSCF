#ifndef _MATH_HPP
#define _MATH_HPP

// Math definitions are placed here
#include <cmath>
#include <Eigen/Dense>
typedef Eigen::Matrix<double, Eigen::Dynamic, 
                              Eigen::Dynamic, Eigen::RowMajor>
    MatrixRowMajor;

// ERI Indexing
#define IJ i*NBF+ j
#define KL k*NBF+ l

#endif // _MATH_HPP