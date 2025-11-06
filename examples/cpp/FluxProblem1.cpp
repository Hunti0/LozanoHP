#include <armadillo>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "utils.h"
#include "gradient.h"
#include "divergence.h"
#include "robinbc.h"
#include "Flux3D.h"

using namespace arma;

// exact solution u(x) = (e^{λ x} - 1) / (e^{λ} - 1)
static vec exact_solution(const vec& x, double lambda) {
    vec u(x.n_elem);
    for (u32 i = 0; i < x.n_elem; ++i)
        u(i) = (std::exp(lambda * x(i)) - 1.0) / (std::exp(lambda) - 1.0);
    return u;
}


static inline double u_x(double x, double lambda) {
    return lambda * std::exp(lambda * x) / (std::exp(lambda) - 1.0);
}

int main() {

    const double lambda = -1.0;
    const double alpha  = -std::exp(lambda);
    const double beta   = (std::exp(lambda) - 1.0) / lambda;

    const u16 k = 2; 
    std::vector<u32> grid_sizes = {50, 100, 200, 400};

    std::vector<double> L2_errors, max_errors;

    std::cout << std::setprecision(10) << std::fixed;
    std::cout << "|   m   |    dx    |   L2 Error   |  L2 Order  |  Max Error  | Max Order |\n";
    std::cout << "|--------|----------|--------------|------------|-------------|-----------|\n";

    for (auto m : grid_sizes) {
        const double dx = 1.0 / static_cast<double>(m);

        
        vec x(m + 2);
        x(0) = 0.0;
        x(m + 1) = 1.0;
        for (u32 i = 1; i <= m; ++i) x(i) = (static_cast<double>(i) - 0.5) * dx;

        vec u_exact = exact_solution(x, lambda);

        // MOLE operators
        Gradient  G(k, m, dx);            
        Divergence Dv(k, m, dx);          
        RobinBC  BC(k, m, dx, alpha, beta); 

        vec kappa(m + 2, fill::ones);
        kappa *= (lambda * lambda);
        auto K = TensorField<1>::Isotropic(kappa);
        FluxND<1> flux(K);

        sp_mat L = flux.diffusion(k, m, dx);

        // Left-hand side and RHS
        sp_mat A = -L + static_cast<sp_mat>(BC);

        vec f = - (sp_mat)L * u_exact;
        f(0)     = alpha * u_exact(0)     - beta * u_x(x(0),     lambda);
        f(m + 1) = alpha * u_exact(m + 1) + beta * u_x(x(m + 1), lambda);

        // Solve
        vec U;
        if (!spsolve(U, A, f, "superlu")) {
            std::cerr << "Solve failed for m = " << m << "\n";
            continue;
        }

        // Errors
        double L2_error  = std::sqrt(dx * sum(square(U - u_exact)));
        double max_error = max(abs(U - u_exact));

        L2_errors.push_back(L2_error);
        max_errors.push_back(max_error);

        std::cout << "| " << std::setw(6) << m
                  << " | " << std::setw(8) << dx
                  << " | " << std::setw(12) << L2_error;

        if (L2_errors.size() > 1) {
            uword i = L2_errors.size() - 1;
            double rate_L2  = std::log(L2_errors[i - 1] / L2_errors[i]) / std::log(2.0);
            double rate_max = std::log(max_errors[i - 1] / max_errors[i]) / std::log(2.0);
            std::cout << " | " << std::setw(10) << rate_L2
                      << " | " << std::setw(11) << max_errors[i]
                      << " | " << std::setw(9)  << rate_max;
        } else {
            std::cout << " | " << std::setw(10) << "-"
                      << " | " << std::setw(11) << max_error
                      << " | " << std::setw(9)  << "-";
        }
        std::cout << " |\n";
    }

    return 0;
}
