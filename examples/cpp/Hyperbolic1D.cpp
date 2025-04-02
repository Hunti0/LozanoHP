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
     
    vec grid(m + 2); //1D Staggered Grid
    grid(0) = west;
    grid(1) = west + dx / 2.0;

    
    for (int i = 2; i <= m; i++) {
        grid(i) = grid(i - 1) + dx;
    }
    grid(m+1) = east;
    
    BC1D periodicBC(k, m, 1, 1, 0); // Create Boundary Conditions

     mat D = -a * dt * 2  * Div * Inter

    return 0;
}
