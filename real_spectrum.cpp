#define _USE_MATH_DEFINES
#include "real_spectrum.h"
#include <cmath>

std::vector<double> realspec::spectrum(double N_droplet, double gdd, double V, double L, int scatter_levels)
{
    int m = 1;
    double E = 0;
    std::vector<double> full_spectrum;

    // First bound level is a droplet level
    full_spectrum.push_back(0);

    // Bound levels
    while (true)
    {
        E = (m / N_droplet) * (9.0 / (4.0 * M_PI * M_PI)) *
            sqrt(0.333333 + 0.25 * (m * m) / (N_droplet * N_droplet)) * gdd * gdd;
        if (E - V < 1e-6)
        {
            full_spectrum.push_back(E);
            m++;
        }
        else
            break;
    }

    // Scatter levels
    for (int k = 1; k <= scatter_levels; k++)
    {
        E = 0.5 * (2 * M_PI * k / L) * (2 * M_PI * k / L) + V;
        full_spectrum.push_back(E);
    }

    return full_spectrum;
}