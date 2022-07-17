#ifndef REAL_SPECTRUM_H
#define REAL_SPECTRUM_H
#include <vector>

class realspec
{
public:
    std::vector<double> spectrum(double N_droplet, double gdd, double V, double L, int scatter_levels);
};

#endif // REAL_SPECTRUM_H