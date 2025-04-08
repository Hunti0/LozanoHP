#include <iostream>
#include <fstream>
#include <cstdlib>
#include <armadillo>
#include "mole.h"

using namespace std;
using namespace arma;

int main() {
    double a = 0.0; // west
    double b = 40.0; // east
    double c = 0.0; // south
    double d = 40.0; // north

    // Grid resolution
    int k = 2; // order of accuracy
    int m = 500; // grid points along x
    int n = 500; // grid points along y
    double dx = (b - a) / m; // Grid spacing in x
    double dy = (d - c) / n; // Grid spacing in y

    // Wave properties
    double wn = 6.0; // Wave number
    double aa = 0.7; // Reflection coefficient
    double bb = 0.9 * 0.5; // Absorption coefficient

    // Hot Spot position
    double hsx = 2.0;
    double hsy = 10.0;
    double hsr = 1.0;

    // 2D Staggered Grid
    vec xgrid = linspace(a, b, m + 2);
    vec ygrid = linspace(c, d, n + 2);
    mat X, Y;

    Util utils;
    utils.mechgrid(xgrid, ygrid, X, Y);

    // Define the Walls
    mat wall = ((X >= 10.0 && X <= 39.0 && Y >= 20.0 && Y <= 21.0) ||
                (X >= 30.0 && X <= 31.0 && Y >= 1.0 && Y <= 16.0) ||
                (X <= 0.5 || X >= 39.5 || Y <= 0.5 || Y >= 39.5));

    // Create Complex Coefficient c = k^2/n^2
    cx_mat ce = pow(wn, 2) / square(1.0 + cx_double(aa, bb) * wall);
    cx_vec ce_vec = vectorise(ce); // Converting to 1D vector to solve

    // Detect hotspot position
    mat HS = (square(X - hsx) + square(Y - hsy) < hsr * hsr);
    vec HS_vec = vectorise(HS); // Convert to 1D vector
    uvec ind = find(HS_vec > 0);
    uvec freenodes = regspace<uvec>(0, (m + 2) * (n + 2) - 1);
    freenodes.shed_rows(ind);

    // Mimetic Laplacian Operator
    Lap2D L(k, m, dx, n, dy);
    L.diag() += ce_vec;

    L += RobinBC2D(k, m, dx, n, dy, 0, 1);

    // RHS
    cx_vec RHS = zeros<cx_vec>((m + 2) * (n + 2));
    RHS -= L * HS_vec; // Modify due to source term

    // Solving Helmholtz
    cx_vec SOL = zeros<cx_vec>((m + 2) * (n + 2));
    SOL(freenodes) = solve(L(freenodes, freenodes), RHS(freenodes));
    SOL(ind).fill(cx_double(1.0, 0.0));

    // Convert back to 2D
    cx_mat SOL2D = reshape(SOL, m + 2, n + 2);

    // Saving solution as log|u|
    mat logSOL = log(abs(SOL2D));
    logSOL.save("solution.dat", raw_ascii); // Save in matrix format

    // Plotting with Gnuplot
    ofstream gp("plot.gp");
    gp << "set title 'Helmholtz Solution (log|u|)'\n";
    gp << "set xlabel 'X'\n";
    gp << "set ylabel 'Y'\n";
    gp << "set view map\n";
    gp << "set pm3d\n";
    gp << "unset key\n";
    gp << "splot 'solution.dat' matrix with image\n";
    gp.close();

    system("gnuplot -persistent plot.gp");

    return 0;
}


