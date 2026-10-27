#include "conjugate_gradient.hpp"

#include <iostream>
#include <iomanip>

int main() {
    Matrix A = {
        {4.0, 1.0},
        {1.0, 3.0}
    };

    Vector b = {1.0, 2.0};
    Vector x0 = {0.0, 0.0};

    CGResult result = conjugateGradient(A, b, x0);

    std::cout << std::fixed << std::setprecision(8);

    std::cout << "Solution:\n";
    for (double value : result.x) {
        std::cout << value << "\n";
    }

    std::cout << "\nIterations: " << result.iterations << "\n";
    std::cout << "Converged: " << std::boolalpha
              << result.converged << "\n";

    std::cout << "\nRelative residual history:\n";
    for (size_t k = 0; k < result.residuals.size(); ++k) {
        std::cout << k << ": " << result.residuals[k] << "\n";
    }

    return 0;
}