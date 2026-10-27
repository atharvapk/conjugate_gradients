
#ifndef CONJUGATE_GRADIENT_HPP
#define CONJUGATE_GRADIENT_HPP

#include <vector>

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct CGResult {
    Vector x;
    std::vector<double> residuals;
    int iterations;
    bool converged;
};

double dotProduct(const Vector& a, const Vector& b);

Vector matrixVectorMultiply(const Matrix& A, const Vector& x);

CGResult conjugateGradient(
    const Matrix& A,
    const Vector& b,
    const Vector& x0,
    double tolerance = 1e-8,
    int maxIterations = 1000
);

#endif