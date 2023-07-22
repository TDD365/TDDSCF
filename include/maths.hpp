#ifndef _MATH_HPP
#define _MATH_HPP

// Math definitions are placed here
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::Matrix<double, Eigen::Dynamic, 
                              Eigen::Dynamic, Eigen::RowMajor>
    MatrixRowMajor;
typedef Eigen::Tensor<double,4,Eigen::RowMajor>
    TensorRowMajor;

#endif // _MATH_HPP