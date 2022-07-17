#define _USE_MATH_DEFINES

#include "Spectrum.h"
#include "well_spectrum.h"
#include "real_spectrum.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

int main()
{
    double gdd = 1000.;                                 // gdd coefficient
    int N_tot = 200;                                    // Number of bosons
    double N_droplet = static_cast<double>(N_tot);      // Number of bosons whose create droplet
    double a = 2.0 * M_PI * M_PI * N_tot / (3.0 * gdd); // Initial droplet width
    double V = 3.0 / (8.0 * M_PI * M_PI) * gdd * gdd;   // Potential
    double L = 10. * a;                                 // Walls
    double c = 2.0 * V * a * a;                         // Constant coefficient in model

    int steps = 20000;                 // Number of sub-bisection steps
    double eps = 1e-6;                 // Accuracy of sub-bisection method
    double xp = eps;                   // Begin of range
    double xkb = sqrt(c * 0.25) - eps; // End of range for bound levels
    double xks = 100.0 - eps;          // End of range for scatter levels

    vector<double> oldSpec;
    vector<double> newSpec;

    Spectrum *S = new Spectrum(eps, steps);
    wellspec *ws = new wellspec(eps, steps);
    realspec *rs = new realspec;

    a = 2.0 * M_PI * M_PI * N_droplet / (3.0 * gdd);
    c = 2.0 * V * a * a;
    xkb = sqrt(c * 0.25) - eps;

    // newSpec = S->calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks);
    // oldSpec = ws->spectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks);

    unsigned scatter_levels = 500; // Number of scatter levels from real spectrum
    oldSpec = rs->spectrum(N_droplet, gdd, V, L, scatter_levels);
    newSpec = S->calculateMasterThesisSpectrum(N_droplet, gdd, V, L, scatter_levels);

    ofstream olds("oldSpec.txt");
    ofstream news("newSpec.txt");

    for (const double &v : newSpec)
        news << v << '\n';

    for (const double &v : oldSpec)
        olds << v << '\n';

    cout << V << '\n';
}