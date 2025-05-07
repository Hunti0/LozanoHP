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
    double dx = (b - a) / m;
    double dy = (d - c) / n;

    // Wave properties
    double wn = 6.0; // wave number
    double aa = 0.7; // reflection coefficient
    double bb = 0.9 * 0.5; // absorption coefficient

    // Hot Spot
    double hsx = 2.0;
    double hsy = 10.0;
    double hsr = 1.0;

    // 2D Staggered Grid
    vec xgrid = linspace(a, b, m + 2);
    vec ygrid = linspace(c, d, n + 2);
    mat X, Y;

    Util utils;
    utils.mechgrid(xgrid, ygrid, X, Y);

    // Define the Walls using element-wise logic
    mat wall = ( ( (X >= 10.0) % (X <= 39.0) % (Y >= 20.0) % (Y <= 21.0) ) +
                 ( (X >= 30.0) % (X <= 31.0) % (Y >= 1.0) % (Y <= 16.0) ) +
                 ( (X <= 0.5) + (X >= 39.5) + (Y <= 0.5) + (Y >= 39.5) ) );
    wall = clamp(wall, 0.0, 1.0); // ensure values are 0 or 1

    // Create complex coefficient c = k^2/n^2
    cx_mat ce = pow(wn, 2) / square(1.0 + cx_double(aa, bb) * wall);
    cx_vec ce_vec = vectorise(ce);

    // Define Hotspot
    mat HS = (square(X - hsx) + square(Y - hsy) < hsr * hsr);
    vec HS_vec_real = vectorise(HS);
    cx_vec HS_vec = conv_to<cx_vec>::from(HS_vec_real); // convert to complex

    // Find hotspot indices
    uvec ind = find(HS_vec_real > 0);
    uvec freenodes = regspace<uvec>(0, (m + 2) * (n + 2) - 1);
    freenodes.shed_rows(ind);

    // Mimetic Laplacian
    Lap2D L(k, m, dx, n, dy);
    L.diag() += ce_vec;

    L += RobinBC2D(k, m, dx, n, dy, 0, 1); // Robin boundary conditions

    // Right-hand side
    cx_vec RHS = zeros<cx_vec>((m + 2) * (n + 2));
    RHS -= L * HS_vec;

    // Solution
    cx_vec SOL = zeros<cx_vec>((m + 2) * (n + 2));
    SOL(freenodes) = solve(L(freenodes, freenodes), RHS(freenodes));
    SOL(ind).fill(cx_double(1.0, 0.0)); // set hotspot manually

    // Reshape solution back to 2D
    cx_mat SOL2D = reshape(SOL, m + 2, n + 2);

    // Save |u| in log scale
    mat logSOL = log(abs(SOL2D));
    logSOL.save("solution.dat", raw_ascii);

    // Gnuplot plotting
    ofstream gp("plot.gp");
    gp << "set title 'Helmholtz Solution'\n";
    gp << "set xlabel 'X'\n";
    gp << "set ylabel 'Y'\n";
    gp << "set size ratio -1\n";
    gp << "unset surface\n";
    gp << "set pm3d map\n";
    gp << "unset key\n";
    gp << "set palette defined ( 0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red' )\n";
    gp << "splot 'solution.dat' matrix with image\n";
    gp.close();

    system("gnuplot -persistent plot.gp");

    return 0;
}
