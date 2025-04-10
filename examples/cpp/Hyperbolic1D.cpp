#include <iostream>
#include <fstream>
#include "mole.h"
#include <vector>
#include <cmath>

using namespace std;

int main() {
    double a = 1; // Velocity
    double west = 0; // Left domain limit
    double east = 1; // Right domain limit

    double k = 2; // Order of Accuracy
    double m = 50; // Number of grid cells

    double dx = (east - west) / m;
    double t = 1; // Simulation Time
    double dt = dx / std::abs(a);

    mat Div = D(k, m, dx);
    mat Inter = I(m, 0.5);

    vec grid(m + 2); // 1D Staggered Grid
    grid(0) = west;
    grid(1) = west + dx / 2.0;
    for (int i = 2; i <= m; i++) {
        grid(i) = grid(i - 1) + dx;
    }
    grid(m + 1) = east;

    // Apply periodic boundary conditions to Div
    Div(0, 1) = 1.0 / (2.0 * dx);
    Div(0, m) = -1.0 / (2.0 * dx);
    Div(m + 1, 1) = 1.0 / (2.0 * dx);
    Div(m + 1, m) = -1.0 / (2.0 * dx);

    mat A = -a * dt * 2 * Div * Inter;

    // Initial condition
    VectorXd U(m + 2);
    for (int i = 0; i < m + 2; ++i) {
        U(i) = sin(2.0 * M_PI * grid(i));
    }

    VectorXd U2 = U + 0.5 * A * U;

    int steps = static_cast<int>(t / dt);

    // Create output file
    ofstream outfile("solution.dat");
    if (!outfile.is_open()) {
        cerr << "Error opening output file!" << endl;
        return 1;
    }

    // Write initial condition
    outfile << "# Time = 0.0" << endl;
    for (int i = 0; i < m + 2; ++i) {
        outfile << grid(i) << " " << U(i) << endl;
    }
    outfile << endl << endl;

    for (int i = 1; i <= steps; ++i) {
        double current_time = i * dt;
        
        // Exact solution for comparison
        VectorXd exact(m + 2);
        for (int j = 0; j < m + 2; ++j) {
            exact(j) = sin(2.0 * M_PI * (grid(j) - a * current_time));
        }

        // Leapfrog time stepping with periodic BC
        VectorXd U3 = U + A * U2;

        // Apply periodic BC after each time step
        U3(0) = U3(m);
        U3(m + 1) = U3(1);

        // Write to file every 10 steps (adjust as needed)
        if (i % 10 == 0) {
            outfile << "# Time = " << current_time << endl;
            for (int j = 0; j < m + 2; ++j) {
                outfile << grid(j) << " " << U3(j) << " " << exact(j) << endl;
            }
            outfile << endl << endl;
        }

        U = U2;
        U2 = U3;
    }

    outfile.close();

    // Generate gnuplot script
    ofstream gnuplot("plot_solution.gp");
    gnuplot << "set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 800, 600\n";
    gnuplot << "set output 'solution.png'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'u(x,t)'\n";
    gnuplot << "set title 'Advection Equation Solution'\n";
    gnuplot << "set key right top\n";
    gnuplot << "plot 'solution.dat' index 0 u 1:2 w l title 'Initial', \\\n";
    gnuplot << "     'solution.dat' index 1 u 1:2 w l title 'Numerical t=0.2', \\\n";
    gnuplot << "     'solution.dat' index 1 u 1:3 w l dt 2 title 'Exact t=0.2', \\\n";
    gnuplot << "     'solution.dat' index 2 u 1:2 w l title 'Numerical t=0.4', \\\n";
    gnuplot << "     'solution.dat' index 2 u 1:3 w l dt 2 title 'Exact t=0.4'\n";
    gnuplot.close();

    return 0;
}
