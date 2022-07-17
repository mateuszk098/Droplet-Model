#define _USE_MATH_DEFINES
#define THREADS 12

#include <complex>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "droplet.h"

using namespace std;

/**
 * @param int total number of bosons
 * @param double initial width of droplet (at T = 0)
 * @param double potential
 **/
Droplet::Droplet(const unsigned &NTot, const double &Width, const double &V) : Ntot(NTot), width(Width), systemPotential(V),
                                                                               integralLowerLimit(-M_PI), integralUpperLimit(M_PI), saddlePointAccuracy(1e-6)
{
    scaleCoefficient = width / Ntot;
    integrationSteps = Ntot * 100;
    integrationAccuracy = (integralUpperLimit - integralLowerLimit) / integrationSteps;
}

/** This function is used to set up energy dropletEnergySpectrum to calculate statistical ensemble.
 * @param std::vector<double> dropletEnergySpectrum (Droplet energy dropletEnergySpectrum)
 * **/
void Droplet::setSpectrum(const std::vector<double> &Spectrum) { dropletEnergySpectrum = Spectrum; }

/** This function calculates grand partition function in grand canonical ensemble
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @return std::complex<double> - gSS (grand partition function)
 * **/
inline complex<double> Droplet::grandStatisticalSum(const double &phi, const double &beta, const double &saddlePoint) const
{
    // We assume that z = exp(u) and z is complex.
    // However, z is related to a physical parameter - the chemical potential u,
    // so the saddle point should exist for a real value of z.
    // I want to integrate along the circle z = r exp(it)
    // which crosses the real saddle point r.

    // Integration trajectory
    complex<double> z = polar(saddlePoint, phi);
    // Grand partition function
    complex<double> gSS = 1.0;

    // Calculation grand partition function for the whole energy dropletEnergySpectrum
    for (const double &energy : dropletEnergySpectrum)
        gSS = gSS * (1.0 / (1.0 - z * exp(-beta * energy)));

    // Return formula for grand partition function
    return gSS;
}

/** This function calculates statistical sum in canonical ensemble based on grand
 * statistical sum using path integral forumaltion.
 * @param std::complex<double> gSS (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @return std::complex<double> - PF (partition function)
 * **/
inline complex<double> Droplet::statisticalSum(const std::complex<double> &gSS, const double &phi, const double &beta,
                                               const double &saddlePoint) const
{
    // Integration trajectory
    complex<double> z = polar(saddlePoint, phi);

    // Return formula for partition function
    return 1.0 / (2.0 * M_PI) * gSS / pow(z, Ntot + 1) * z;
}

/** This function calculates average stateOccupation of a specific state in canonical ensemble
 * based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> gSS (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @param int level (energy dropletEnergySpectrum level)
 * @return std::complex<double> - AO (average stateOccupation of a specific state)
 * **/
inline complex<double> Droplet::stateOccupation(const std::complex<double> &gSS, const double &phi, const double &beta,
                                                const double &saddlePoint, const unsigned &specState) const
{
    // Integration trajectory
    complex<double> z = polar(saddlePoint, phi);
    // Energy of a specific state
    double specific_level = dropletEnergySpectrum[specState];

    // Return formula for stateOccupation
    return 1.0 / (2.0 * M_PI) * z * exp(-beta * specific_level) / (1.0 - z * exp(-beta * specific_level)) * gSS / pow(z, Ntot + 1) * z;
}

/** This function calculates stateFluctuations of average stateOccupation of a specific state in canonical
 * ensemble based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> gSS (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double saddle_point (saddle point for the integration process)
 * @param int level (energy dropletEnergySpectrum level)
 * @return std::complex<double> - FL (stateFluctuations of a specific state)
 * **/
inline complex<double> Droplet::stateFluctuations(const std::complex<double> &gSS, const double &phi, const double &beta,
                                                  const double &saddlePoint, const unsigned &specState) const
{
    // Integration trajectory
    complex<double> z = polar(saddlePoint, phi);
    // Energy of a specific state
    double specific_level = dropletEnergySpectrum[specState];

    // Return formula for stateFluctuations
    return 1.0 / (2 * M_PI) * z * exp(-beta * specific_level) / (1.0 - z * exp(-beta * specific_level)) * z * exp(-beta * specific_level) / (1.0 - z * exp(-beta * specific_level)) * gSS / pow(z, Ntot + 1) * z;
}

/** This function is used to calculate function to saddle point.
 * @param double z (searched saddle point)
 * @param double beta (beta coefficient related to temperature)
 * @param int Ntot (number of bosons)
 * @return double - function to saddle point
 * **/
inline double Droplet::saddlePointFunction(const double &z, const double &beta, const unsigned &Ntot) const
{
    // Function on a saddle point
    double func = 0.;

    for (const double &energy : dropletEnergySpectrum)
        func += z / (exp(beta * energy) - z);

    return func - (Ntot + 1);
}

/** Bisection method to find saddle point to integration process
 * Unfortunately it does not work as well as for condensate so we are unlikely to use it.
 * Perhaps we will be able to use this but within in narrow range of a and b, e.g. a = 0.9
 * and b = 1.0. In any case we have to use Mathematica for low temperatures, where saddle
 * point is higher than 1.0
 * @param double start (start point of a range)
 * @param double end (end point of a range)
 * @param double beta (beta coefficient (beta = 1/T))
 * @param int Ntot (total number of bosons)
 * @return double z0 (saddle point for the integration process)
 * **/
double Droplet::saddlePoint(double start, double end, const double &beta, const unsigned &Ntot) const
{
    double z0 = 0.; // Searched saddle point for the integration process

    if (saddlePointFunction(start, beta, Ntot) * saddlePointFunction(end, beta, Ntot) > 0.)
        return -1.;
    else
    {
        while ((end - start) * 0.5 > saddlePointAccuracy)
        {
            z0 = (start + end) * 0.5;

            if (saddlePointFunction(z0, beta, Ntot) == 0)
                break;
            else if (saddlePointFunction(start, beta, Ntot) * saddlePointFunction(z0, beta, Ntot) < 0.)
                end = z0;
            else
                start = z0;
        }
    }

    return z0;
}

/** This function is used to calculate the specific state stateOccupation and
 * its basic properties. In this function we set beta coefficient and energy dropletEnergySpectrum level.
 * @param double beta (beta coefficient (beta = 1/T))
 * @param int specState (specific state of energy dropletEnergySpectrum)
 * @return textfile (text file with properties)
 * **/
void Droplet::stateProperties(const double &T, const unsigned &specState, std::string filename)
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile.precision(5);

    double beta = 1 / T;
    double saddle_point = saddlePoint(0.0, 1.0, beta, Ntot); // Saddle point for the integration process
    double phi = 0;                                          // Integration step
    complex<double> gSS = 0;                                 // Grand partition function
    double PF = 0;                                           // Partition function
    double AO = 0;                                           // Average stateOccupation
    double FL = 0;                                           // Fluctuations

#pragma omp parallel for private(phi, gSS) reduction(+ \
                                                     : PF, AO, FL) num_threads(THREADS)
    for (int k = 1; k <= integrationSteps; k++)
    {
        phi = integralLowerLimit + k * integrationAccuracy;
        gSS = grandStatisticalSum(phi, beta, saddle_point);
        PF = PF + statisticalSum(gSS, phi, beta, saddle_point).real();
        AO = AO + stateOccupation(gSS, phi, beta, saddle_point, specState).real();
        FL = FL + stateFluctuations(gSS, phi, beta, saddle_point, specState).real();
    }

    PF = PF * integrationAccuracy;
    AO = AO / PF * integrationAccuracy;
    FL = 2.0 / PF * FL * integrationAccuracy + AO - AO * AO;

    ofile << "Energy level:                          " << specState << "\n";
    ofile << "Total number of bosons:                " << Ntot << "\n";
    ofile << "Temperature:                           " << 1.0 / beta << "\n";
    ofile << "Saddle point:                          " << saddle_point << "\n";
    ofile << "Partition function:                    " << PF << "\n";
    ofile << "Helmholtz free energy:                 " << -log(PF) << "\n";
    ofile << "Average stateOccupation:                    " << AO << "\n";
    ofile << "Fluctuations of average stateOccupation:    " << FL << "\n";
    ofile << "Dispersion of average stateOccupation:      " << sqrt(FL) << "\n";
}

/** This function is used to check what happen with specific state
 * stateOccupation and stateFluctuations with temperature change.
 * @param int specState (specific state of energy dropletEnergySpectrum)
 * @return Text file with Droplet properties
 * **/
void Droplet::stateTemperatureImpact(const unsigned &specState, const double &Tp, const double &Tk, const double &dT, std::string filename)
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile.precision(5);
    ofile << std::fixed;
    ofile << std::showpos;

    double beta = 0;         // Beta coefficient
    double saddle_point = 0; // Real saddle point
    double phi = 0;          // Integration step

    complex<double> gSS = 0; // Grand partition function
    double PF = 0;           // Partition function
    double AO = 0;           // Average stateOccupation
    double FL = 0;           // Fluctuations

    ofile << "T\t"
          << "F\t"
          << "<N>\t"
          << "Sigma\t"
          << "\n";

    for (double T = Tp; T <= Tk; T += dT)
    {
        beta = 1.0 / T;                                   // Beta coefficient
        saddle_point = saddlePoint(0.0, 1.0, beta, Ntot); // Real saddle point

#pragma omp parallel for private(phi, gSS) reduction(+ \
                                                     : PF, AO, FL) num_threads(THREADS)
        for (int k = 1; k <= integrationSteps; k++)
        {
            phi = integralLowerLimit + k * integrationAccuracy;
            gSS = grandStatisticalSum(phi, beta, saddle_point);
            PF = PF + statisticalSum(gSS, phi, beta, saddle_point).real();
            AO = AO + stateOccupation(gSS, phi, beta, saddle_point, specState).real();
            FL = FL + stateFluctuations(gSS, phi, beta, saddle_point, specState).real();
        }

        PF = PF * integrationAccuracy;
        AO = AO / PF * integrationAccuracy;
        FL = 2.0 / PF * FL * integrationAccuracy + AO - AO * AO;

        ofile << T << "\t" << -log(PF) << "\t" << AO << "\t" << FL << "\n";
    }

    ofile.close();
}

/** This function is use to calculate Droplet width in a special temperature.
 * Quantum Droplet consist of first Droplet state and every bound states
 * @param double beta (beta coefficient (beta = 1/T))
 * @param double& N_droplet (number of bosons whose create Droplet)
 * @param bool flag (flag of text file save)
 * @return Text file with Droplet properties
 * **/
std::stringstream Droplet::dropletWidth(const double &T, double &NDroplet)
{
    NDroplet = 0; // Number of bosons whose create quantum Droplet
    double beta = 1.0 / T;
    double saddle_point = saddlePoint(0.0, 1.0, beta, Ntot);
    double Fluctuations = 0; // Fluctuations of every bound states
    double F = 0;            // Helmholtz free energy
    int max_bound_level = 0; // Number of every states whose energy is less than potential

    // Counting bound levels
    while (dropletEnergySpectrum[max_bound_level] - systemPotential < saddlePointAccuracy)
        max_bound_level++;

    double phi;              // Integration step
    complex<double> gSS = 0; // Grand partition function
    double PF = 0;           // Partition function
    double AO = 0;           // Average stateOccupation
    double FL = 0;           // Fluctuations

    for (int specState = 0; specState < max_bound_level; specState++)
    {
#pragma omp parallel for private(phi, gSS) reduction(+ \
                                                     : PF, AO, FL) num_threads(THREADS)
        for (int k = 1; k <= integrationSteps; k++)
        {
            phi = integralLowerLimit + k * integrationAccuracy;
            gSS = grandStatisticalSum(phi, beta, saddle_point);
            PF = PF + statisticalSum(gSS, phi, beta, saddle_point).real();
            AO = AO + stateOccupation(gSS, phi, beta, saddle_point, specState).real();
            FL = FL + stateFluctuations(gSS, phi, beta, saddle_point, specState).real();
        }

        PF = PF * integrationAccuracy;
        AO = AO / PF * integrationAccuracy;
        FL = 2.0 / PF * FL * integrationAccuracy + AO - AO * AO;

        NDroplet += AO;
        Fluctuations += FL;
        F += -log(PF);
    }

    stringstream output;
    output << T << "\t" << scaleCoefficient * NDroplet << "\t" << scaleCoefficient * Fluctuations << "\t" << F << "\n";
    return output;
}