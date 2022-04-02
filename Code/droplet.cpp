#define _USE_MATH_DEFINES
#define THREADS 12

#include <complex>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "droplet.h"

using namespace std;

/** This is constructor. It sets up proper values to initial droplet and
 * integration process. Also it calculates and sets alpha coefficient.
 * @param int N_tot (total number of bosons)
 * @param double width (initial width of droplet (at T = 0))
 * @param double V0 (potential (height of initial droplet))
 * **/
droplet::droplet(int this_N_tot, double this_width, double this_V0)
{
    // Initial droplet parameters
    N_tot = this_N_tot; // Number of bosons
    V0 = this_V0;       // Number ofPotential
    width = this_width; // Droplet width
    alpha = width / N_tot;

    // Path intergration parameters
    tp = -M_PI;         // Lower limit of the integral
    tk = M_PI;          // Upper limit of the integral
    n = N_tot * 100;    // Number of interation steps should be 100-time greater
    dt = (tk - tp) / n; // Integration step

    // Saddle point parameters
    eps = 1e-6; // Saddle point accuracy
}

/** This function is used to set up energy spectrum to calculate statistical ensemble.
 * @param std::vector<double> spectrum (droplet energy spectrum)
 * **/
void droplet::set_spectrum(vector<double> this_spectrum)
{
    spectrum = this_spectrum; // Energy spectrum, first state is a droplet state
}

/** This function calculates grand partition function in grand canonical ensemble
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @return std::complex<double> - GPF (grand partition function)
 * **/
complex<double> droplet::grandstatsum(double phi, double beta, double saddle_point)
{
    // We assume that z = exp(u) and z is complex.
    // However, z is related to a physical parameter - the chemical potential u,
    // so the saddle point should exist for a real value of z.
    // I want to integrate along the circle z = r exp(it)
    // which crosses the real saddle point r.

    // Integration trajectory
    complex<double> z = polar(saddle_point, phi);
    // Grand partition function
    complex<double> GPF = 1.0;

    // Calculation grand partition function for the whole energy spectrum
    for (auto &energy : spectrum)
        GPF = GPF * (1.0 / (1.0 - z * exp(-beta * energy)));

    // Return formula for grand partition function
    return GPF;
}

/** This function calculates statistical sum in canonical ensemble based on grand
 * statistical sum using path integral forumaltion.
 * @param std::complex<double> GPF (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @return std::complex<double> - PF (partition function)
 * **/
complex<double> droplet::statsum(complex<double> GPF, double phi, double beta, double saddle_point)
{
    // Integration trajectory
    complex<double> z = polar(saddle_point, phi);

    // Return formula for partition function
    return 1.0 / (2.0 * M_PI) * GPF / pow(z, N_tot + 1) * z;
}

/** This function calculates average occupation of a specific state in canonical ensemble
 * based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> GPF (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @param int level (energy spectrum level)
 * @return std::complex<double> - AO (average occupation of a specific state)
 * **/
complex<double> droplet::occupation(complex<double> GPF, double phi, double beta, double saddle_point, int level)
{
    // Integration trajectory
    complex<double> z = polar(saddle_point, phi);
    // Energy of a specific state
    double specific_level = spectrum[level];

    // Return formula for occupation
    return 1.0 / (2.0 * M_PI) * z * exp(-beta * specific_level) / (1.0 - z * exp(-beta * specific_level)) * GPF / pow(z, N_tot + 1) * z;
}

/** This function calculates fluctuations of average occupation of a specific state in canonical
 * ensemble based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> GPF (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @param int level (energy spectrum level)
 * @return std::complex<double> - FL (fluctuations of a specific state)
 * **/
complex<double> droplet::fluctuations(complex<double> GPF, double phi, double beta, double saddle_point, int level)
{
    // Integration trajectory
    complex<double> z = polar(saddle_point, phi);
    // Energy of a specific state
    double specific_level = spectrum[level];

    // Return formula for fluctuations
    return 1.0 / (2 * M_PI) * z * exp(-beta * specific_level) / (1.0 - z * exp(-beta * specific_level)) * z * exp(-beta * specific_level) / (1.0 - z * exp(-beta * specific_level)) * GPF / pow(z, N_tot + 1) * z;
}

/** This function is used to calculate function to saddle point.
 * @param double z (searched saddle point)
 * @param double beta (beta coefficient related to temperature)
 * @param int N_tot (number of bosons)
 * @return double - function to saddle point
 * **/
double droplet::SPE(double z, double beta, int N_tot)
{
    // Function on a saddle point
    double func = 0;

    for (auto &energy : spectrum)
        func += z / (exp(beta * energy) - z);

    return func - (N_tot + 1);
}

/** Bisection method to find saddle point to integration process
 * Unfortunately it does not work as well as for condensate so we are unlikely to use it.
 * Perhaps we will be able to use this but within in narrow range of a and b, e.g. a = 0.9
 * and b = 1.0. In any case we have to use Mathematica for low temperatures, where saddle
 * point is higher than 1.0
 * @param double start (start point of a range)
 * @param double end (end point of a range)
 * @param double beta (beta coefficient (beta = 1/T))
 * @param int N_tot (total number of bosons)
 * @return double z0 (saddle point for the integration process)
 * **/
double droplet::SaddlePoint(double start, double end, double beta, int N_tot)
{
    double z0 = 0; // Searched saddle point for the integration process

    if (SPE(start, beta, N_tot) * SPE(end, beta, N_tot) > 0)
        return -1;
    else
    {
        while ((end - start) * 0.5 > eps)
        {
            z0 = (start + end) * 0.5;

            if (SPE(z0, beta, N_tot) == 0)
                break;
            else if (SPE(start, beta, N_tot) * SPE(z0, beta, N_tot) < 0)
                end = z0;
            else
                start = z0;
        }
    }

    return z0;
}

/** This function is used to calculate the specific state occupation and
 * its basic properties. In this function we set beta coefficient and energy spectrum level.
 * @param double beta (beta coefficient (beta = 1/T))
 * @param int spectrum_level (specific state of energy spectrum)
 * @return textfile (text file with properties)
 * **/
void droplet::specific_state_properties(double T, int spectrum_level, string filename)
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile.precision(5);

    double beta = 1 / T;
    double saddle_point = SaddlePoint(0.0, 1.0, beta, N_tot); // Saddle point for the integration process
    double phi = 0;                                           // Integration step
    complex<double> GPF = 0;                                  // Grand partition function
    double PF = 0;                                            // Partition function
    double AO = 0;                                            // Average occupation
    double FL = 0;                                            // Fluctuations

#pragma omp parallel for private(phi, GPF) reduction(+ \
                                                     : PF, AO, FL) num_threads(THREADS)
    for (int k = 1; k <= n; k++)
    {
        phi = tp + k * dt;
        GPF = grandstatsum(phi, beta, saddle_point);
        PF = PF + statsum(GPF, phi, beta, saddle_point).real();
        AO = AO + occupation(GPF, phi, beta, saddle_point, spectrum_level).real();
        FL = FL + fluctuations(GPF, phi, beta, saddle_point, spectrum_level).real();
    }

    PF = PF * dt;
    AO = AO / PF * dt;
    FL = 2.0 / PF * FL * dt + AO - AO * AO;

    ofile << "Energy level:                          " << spectrum_level << "\n";
    ofile << "Total number of bosons:                " << N_tot << "\n";
    ofile << "Temperature:                           " << 1.0 / beta << "\n";
    ofile << "Saddle point:                          " << saddle_point << "\n";
    ofile << "Partition function:                    " << PF << "\n";
    ofile << "Helmholtz free energy:                 " << -log(PF) << "\n";
    ofile << "Average occupation:                    " << AO << "\n";
    ofile << "Fluctuations of average occupation:    " << FL << "\n";
    ofile << "Dispersion of average occupation:      " << sqrt(FL) << "\n";
}

/** This function is used to check what happen with specific state
 * occupation and fluctuations with temperature change.
 * @param int spectrum_level (specific state of energy spectrum)
 * @return Text file with droplet properties
 * **/
void droplet::state_with_temperature_change(int spectrum_level, double Tp, double Tk, double dT, string filename)
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile.precision(5);
    ofile << std::fixed;
    ofile << std::showpos;

    double beta = 0;         // Beta coefficient
    double saddle_point = 0; // Real saddle point
    double phi = 0;          // Integration step

    complex<double> GPF = 0; // Grand partition function
    double PF = 0;           // Partition function
    double AO = 0;           // Average occupation
    double FL = 0;           // Fluctuations

    ofile << "T\t"
          << "F\t"
          << "<N>\t"
          << "Sigma\t"
          << "\n";

    for (double T = Tp; T <= Tk; T += dT)
    {
        beta = 1.0 / T;                                    // Beta coefficient
        saddle_point = SaddlePoint(0.0, 1.0, beta, N_tot); // Real saddle point

#pragma omp parallel for private(phi, GPF) reduction(+ \
                                                     : PF, AO, FL) num_threads(THREADS)
        for (int k = 1; k <= n; k++)
        {
            phi = tp + k * dt;
            GPF = grandstatsum(phi, beta, saddle_point);
            PF = PF + statsum(GPF, phi, beta, saddle_point).real();
            AO = AO + occupation(GPF, phi, beta, saddle_point, spectrum_level).real();
            FL = FL + fluctuations(GPF, phi, beta, saddle_point, spectrum_level).real();
        }

        PF = PF * dt;
        AO = AO / PF * dt;
        FL = 2.0 / PF * FL * dt + AO - AO * AO;

        ofile << T << "\t" << -log(PF) << "\t" << AO << "\t" << FL << "\n";
    }

    ofile.close();
}

/** This function is use to calculate droplet width in a special temperature.
 * Quantum droplet consist of first droplet state and every bound states
 * @param double beta (beta coefficient (beta = 1/T))
 * @param double& N_droplet (number of bosons whose create droplet)
 * @param bool flag (flag of text file save)
 * @return Text file with droplet properties
 * **/
std::stringstream droplet::droplet_width(double T, double &N_droplet)
{
    N_droplet = 0; // Number of bosons whose create quantum droplet
    double beta = 1.0 / T;
    double saddle_point = SaddlePoint(0.0, 1.0, beta, N_tot);
    double Fluctuations = 0; // Fluctuations of every bound states
    double F = 0;            // Helmholtz free energy
    int max_bound_level = 0; // Number of every states whose energy is less than potential

    // Counting bound levels
    while (spectrum[max_bound_level] - V0 < eps)
        max_bound_level++;

    double phi;              // Integration step
    complex<double> GPF = 0; // Grand partition function
    double PF = 0;           // Partition function
    double AO = 0;           // Average occupation
    double FL = 0;           // Fluctuations

    for (int spectrum_level = 0; spectrum_level < max_bound_level; spectrum_level++)
    {
#pragma omp parallel for private(phi, GPF) reduction(+ \
                                                     : PF, AO, FL) num_threads(THREADS)
        for (int k = 1; k <= n; k++)
        {
            phi = tp + k * dt;
            GPF = grandstatsum(phi, beta, saddle_point);
            PF = PF + statsum(GPF, phi, beta, saddle_point).real();
            AO = AO + occupation(GPF, phi, beta, saddle_point, spectrum_level).real();
            FL = FL + fluctuations(GPF, phi, beta, saddle_point, spectrum_level).real();
        }

        PF = PF * dt;
        AO = AO / PF * dt;
        FL = 2.0 / PF * FL * dt + AO - AO * AO;

        N_droplet += AO;
        Fluctuations += FL;
        F += -log(PF);
    }

    stringstream output;
    output << T << "\t" << alpha * N_droplet << "\t" << alpha * Fluctuations << "\t" << F << "\n";
    return output;
}