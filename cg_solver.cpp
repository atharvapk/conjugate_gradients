
#include "conjugate_gradient.hpp"

#include <cmath>
#include <stdexcept>

double dotProduct(const Vector& a, const Vector& b) {
    // Vector dot products only work when the dimensions are a 1:1 match
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector dot products only work when the dimensions are a 1:1 match");
    }

    double result = 0.0;

    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }

    return result;
}

Vector matrixVectorMultiply(const Matrix& A, const Vector& x) {
    if (A.size() != x.size()) {
        throw std::invalid_argument("Matrix and vector dimensions do not match");
    }

    Vector result(A.size(), 0.0);

    for (size_t i = 0; i < A.size(); ++i) {
        // The length of each row in A must match the height of x
        // In this case, this coincides with A needing to be a square
        // This is because if the lenght of a row of A does not match the height of x, A is mxn
        if (A[i].size() != x.size()) {
            throw std::invalid_argument("Matrix must be square");
        }

        for (size_t j = 0; j < x.size(); ++j) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}

CGResult conjugateGradient(
    const Matrix& A,
    const Vector& b,
    const Vector& x0,
    double tolerance,
    int maxIterations
) {
    const size_t n = b.size();

    if (x0.size() != n || A.size() != n) {
        throw std::invalid_argument("Dimension mismatch");
    }
    if (tolerance <= 0.0 || maxIterations < 0) {
        throw std::invalid_argument("Invalid solver parameters");
    }

    Vector x = x0;
    Vector Ax = matrixVectorMultiply(A, x);
    Vector r(n), p(n);

    for (size_t i = 0; i < n; ++i) {
        r[i] = b[i] - Ax[i];
    }

    p = r;

    double rr = dotProduct(r, r);
    double normB = std::sqrt(dotProduct(b, b));
    double scale = (normB > 0.0) ? normB : 1.0;

    std::vector<double> residuals;
    residuals.push_back(std::sqrt(rr) / scale);

    if (residuals.back() <= tolerance) {
        return {x, residuals, 0, true};
    }

    for (int k = 0; k < maxIterations; ++k) {
        Vector Ap = matrixVectorMultiply(A, p);

        double pAp = dotProduct(p, Ap);

        if (pAp <= 0.0) {
            throw std::runtime_error(
                "CG breakdown: matrix may not be positive definite"
            );
        }

        double alpha = rr / pAp;

        for (size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rrNew = dotProduct(r, r);
        double relativeResidual = std::sqrt(rrNew) / scale;

        residuals.push_back(relativeResidual);

        if (relativeResidual <= tolerance) {
            return {x, residuals, k + 1, true};
        }

        double beta = rrNew / rr;

        for (size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        rr = rrNew;
    }

    return {x, residuals, maxIterations, false};
}