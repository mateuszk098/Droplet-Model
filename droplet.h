#ifndef DROPLET_H
#define DROPLET_H

#include <vector>
#include <complex>
#include <sstream>

class droplet
{
private:
    int n;                        // Integration steps
    double tp;                    // Lower limit of the integral
    double tk;                    // Upper limit of the integral
    double dt;                    // Integration step
    int N_tot;                    // Total number of bosons
    double V0;                    // Potential
    double width;                 // Initial droplet width
    double alpha;                 // Constant scale coefficient beetwen N_tot and width
    std::vector<double> spectrum; // Energy spectrum, first state is a droplet state
    double eps;                   // Saddle point accuracy

    std::complex<double> grandstatsum(double phi, double beta, double saddle_point);
    std::complex<double> statsum(std::complex<double> GPF, double phi, double beta, double saddle_point);
    std::complex<double> occupation(std::complex<double> GPF, double phi, double beta, double saddle_point, int spectrum_level);
    std::complex<double> fluctuations(std::complex<double> GPF, double phi, double beta, double saddle_point, int spectrum_level);

    double SPE(double z, double beta, int N_tot);
    double SaddlePoint(double start, double end, double beta, int N_tot);

public:
    droplet(int this_N_tot, double this_width, double this_V0);
    void set_spectrum(std::vector<double> this_spectrum);
    void specific_state_properties(double beta, int spectrum_level);
    void state_with_temperature_change(int spectrum_level);
    std::stringstream droplet_width(double beta, double &N_droplet);
};

#endif // DROPLET_H