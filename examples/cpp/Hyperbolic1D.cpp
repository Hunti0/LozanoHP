#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "mole.h"

using namespace std;
using namespace Eigen;

int main() {
    double a = 1.0;       // Advection velocity
    double west = 0.0;    // Left boundary
    double east = 1.0;    // Right boundary
    double k = 2;         // Order of accuracy
    int m = 50;           // Number of cells

    double dx = (east - west) / m;
    double t = 1.0;       // Final simulation time
    double dt = dx / abs(a); // CFL time step
    int steps = static_cast<int>(t / dt);

    Divergence  D(k, m, dx);
    Interpol  I(m, 0.5);

    // Grid definition: staggered grid with ghost points
    vec grid(m + 2);
    grid(0) = west;
    grid(1) = west + dx / 2.0;
    for (int i = 2; i <= m; ++i) {
        grid(i) = grid(i - 1) + dx;
    }
    grid(m + 1) = east;

    // Apply periodic boundary conditions to Div
    D(0, 1) = 1.0 / (2.0 * dx);
    D(0, m) = -1.0 / (2.0 * dx);
    D(m + 1, 1) = 1.0 / (2.0 * dx);
    D(m + 1, m) = -1.0 / (2.0 * dx);

    mat A = -a * dt * 2.0 * Div * Inter;

    // Initial condition
    VectorXd U(m + 2);
    for (int i = 0; i < m + 2; ++i) {
        U(i) = sin(2.0 * M_PI * grid(i));
    }

    // Bootstrap leapfrog using Euler step
    VectorXd U2 = U + 0.5 * A * U;

    // Open output file
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

    // Time stepping using Leapfrog
    for (int i = 1; i <= steps; ++i) {
        VectorXd U3 = U + A * U2;

        // Periodic BCs
        U3(0) = U3(m);
        U3(m + 1) = U3(1);

        // Output every 10 steps
        if (i % 10 == 0) {
            double current_time = i * dt;
            outfile << "# Time = " << current_time << endl;
            for (int j = 0; j < m + 2; ++j) {
                double x_shifted = grid(j) - a * current_time;

                // Periodicity for exact solution
                while (x_shifted < west) x_shifted += (east - west);
                while (x_shifted > east) x_shifted -= (east - west);

                double exact_val = sin(2.0 * M_PI * x_shifted);
                outfile << grid(j) << " " << U3(j) << " " << exact_val << endl;
            }
            outfile << endl << endl;
        }

        U = U2;
        U2 = U3;
    }

    outfile.close();

    // Gnuplot script
    ofstream gnuplot("plot_solution.gp");
    gnuplot << "set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 800, 600\n";
    gnuplot << "set output 'solution.png'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'u(x,t)'\n";
    gnuplot << "set title 'Advection Equation '\n";
    gnuplot << "set key right top\n";
    gnuplot << "plot 'solution.dat' index 0 u 1:2 w l title 'Initial', \\\n";
    gnuplot << "     'solution.dat' index 1 u 1:2 w l title 'Numerical t=0.2', \\\n";
    gnuplot << "     'solution.dat' index 1 u 1:3 w l dt 2 title 'Exact t=0.2', \\\n";
    gnuplot << "     'solution.dat' index 2 u 1:2 w l title 'Numerical t=0.4', \\\n";
    gnuplot << "     'solution.dat' index 2 u 1:3 w l dt 2 title 'Exact t=0.4'\n";
    gnuplot.close();

    return 0;
}
