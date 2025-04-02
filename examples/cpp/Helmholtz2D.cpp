#include <iostream>
#include <vector>
#include <cmath>
#include "mole.h" 
int main() {
    
const double aa = 0.7;   // Wall reflection coefficient
const double bb = 0.9 * 0.5;  // Wall absorption coefficient
const double wn = 6.0;   // Angular wave number

 //Spatial Discretization
const int k = 2;  // Order of accuracy
const int m = 500, n = 500;  // Grid resolution
const double a = 0, b = 40;  // X domain
const double c = 0, d = 40;  // Y domain
const double dx = (b - a) / m; // Grid spacing
const double dy = (d - c) / n;  // Grid spacing

// Hotspot location
const double hsx = 2.0;
const double hsy = 10.0;
const double hsr = 1.0;

bool isWall(double x, double y) {
    return ((x >= 10.0 && x <= 39.0 && y >= 20.0 && y <= 21.0) ||
            (x >= 30.0 && x <= 31.0 && y >= 1.0 && y <= 16.0) ||
            (x <= 0.5 || x >= 39.5 || y <= 0.5 || y >= 39.5));
}

bool isHotspot(double x, double y) {
    return ((x - hsx) * (x - hsx) + (y - hsy) * (y - hsy) < hsr * hsr);
}


  Lap2D L(k, grid);
   


    return 0;
}
