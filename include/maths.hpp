#ifndef _MATH_HPP
#define _MATH_HPP

// Math definitions are placed here
#include <cmath>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
typedef Eigen::Matrix<double, Eigen::Dynamic, 
                              Eigen::Dynamic, Eigen::RowMajor>
    MatrixRowMajor;
typedef Eigen::Tensor<double,4,Eigen::RowMajor>
    TensorRowMajor;

// ERI Indexing
#define IJ i*NBF+ j
#define KL k*NBF+ l
#define IL i*NBF+ l
#define KJ k*NBF+ j

#endif // _MATH_HPP