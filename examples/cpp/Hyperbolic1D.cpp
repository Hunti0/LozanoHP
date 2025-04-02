#include <iostream>
#include "mole.h"
#include <vector>
#include <cmath>
int main() {
    double a = 1; //Velocity
    double west = 0; // Left domain limit
    double east  = 1; // Right domain limit

    double k = 2; //Order of Accuracy
    double m  = 50; // Number of grid cells

    double dx = (east - west) / m; // Time Step
    int t = 1; //Simulation Time
    double dt = dx / std::abs(a);

    Div = D(k, m, dx);
    Inter = I(m, 0.5);
     //grid 
    BC1D periodicBC(k, m, 1, 1, 0);

     mat D = -a * dt * 2  * Div * Inter

     vector<double> U(m), U2(m), U3(m);
    for (int i = 0; i < m; ++i) {
        double x = grid.x(i);
        U[i] = sin(2 * M_PI * x); // Initial condition
    }

    U2 = U + 0.5 * (Div * U)

    for (double t = dt; t<= t_final; t += dt) {
        U3 = U + (D * U2);
        U = U2;
        U2 = U3;
        
    return 0;
}
