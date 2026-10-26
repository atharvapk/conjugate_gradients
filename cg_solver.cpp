#include <vector>
#include <cmath>
#include <stdexcept>

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct CGResult {
    Vector x;
    std::vector<double> residuals;
    int iterations;
    bool converged;
};

CGResult conjugateGradient(
    const Matrix& A,
    const Vector& b,
    const Vector& x0,
    double tol = 1e-8,
    int maxIter = 1000
) {
    Vector x = x0;
    Vector r(b.size());
    Vector p(b.size());

    //compute r = b-A*X
};