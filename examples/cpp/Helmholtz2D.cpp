#include <iostream>
#include <armadillo>
#include "mole.h"

using namespace std:
using namepace arma;

int main() {
    double a = 0.0; // west
    double b = 40.0; // east
    double c = 0.0; // south
    double d = 40.0; // north
// Grid resolution
int k = 2; //order of accuracy
int m = 500; // grid points along x
int n = 500 // grid points along y
double dx = (b - a) / m; //Grid spacing in x
double dy = (d - c) / n //Grid spacing in y

// Wave properties
double wn = 6.0; // Wave number
double aa = 0.7; // Reflection coefficient
double bb = 0.9 * 0.5; // Absorbtion coefficient

// Hot Spot position
double hsx = 2.0; 
double hsy = 10.0;
double hsr = 1.0;

// Define the Walls
mat wall = ( (X >= 10.0 && X <= 39.0 && Y >= 20.0 && Y <= 21.0) ||
                 (X >= 30.0 && X <= 31.0 && Y >= 1.0 && Y <= 16.0) ||
                 (X <= 0.5 || X >= 39.5 || Y <= 0.5 || Y >= 39.5) );

// 2D Staggered Grid
vec xgrid = linspace(a, b, m + 2);
vec ygrid = linspace(c, d, n + 2);

mat X, Y;

Util utils;
utils.mechgrid(xgrid, ygrid, X, Y);

// Create Complex Coefficient c = k^2/n^2
mat ce = (wn * wn) / pow(1+ (aa + cx_double(0, bb)) * wall, 2);
vec ce_vec = vectorise(ce) // Converting to 1D vector to solve

// Detect hostpot position
mat HS = (square(X- hsx) + square(Y - hsy) < hsr * hsr);
vec HS_vec = vectorise(HS); // Convert to 1D vector
uvec ind = find(HS_vec > 0);
uvec freenodes =  regenspace<uvec>(0, (m+2) * (n+2) - 1);
freenodes.shed_rows(ind);

// Mimetic Laplacian Operator
Laplacian2D L (k, m , dx, n, dy);


