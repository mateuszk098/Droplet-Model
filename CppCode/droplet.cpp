#define _USE_MATH_DEFINES
#define THREADS 12

#include <complex>
#include <cmath>
#include <fstream>
#include <omp.h>
#include "droplet.h"
#include <iomanip>
using namespace std;

Droplet::Droplet(const unsigned &NTot, const double &Width, const double &V) : totalParticlesNumber(NTot), width(Width), systemPotential(V),
                                                                               integralLowerLimit(-M_PI), integralUpperLimit(M_PI), saddlePointAccuracy(1e-6)
{
    scaleCoefficient = width / totalParticlesNumber;
    integrationSteps = totalParticlesNumber * 100;
    integrationAccuracy = (integralUpperLimit - integralLowerLimit) / integrationSteps;
}

void Droplet::setSpectrum(const std::vector<double> &dropletEnergySpectrum) { this->dropletEnergySpectrum = dropletEnergySpectrum; }

/** Calculates grand partition function in grand canonical ensemble
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double foundedSaddlePoint (saddle point for the integration process)
 * @return std::complex<double> - grandStatSum (grand statistical sum)
 * **/
inline complex<double> Droplet::grandStatisticalSum(const double &phi, const double &beta, const double &saddlePoint) const
{
    // We assume that z = exp(u) and z is complex. However, z is related to a physical parameter - the chemical
    // potential u, so the saddle point should exist for a real value of z. I want to integrate along the circle
    // z = r exp(it) which crosses the real saddle point r.

    const complex<double> z = polar(saddlePoint, phi);
    complex<double> grandStatSum = 1.;

    for (const double &energy : dropletEnergySpectrum)
        grandStatSum = grandStatSum * (1. / (1. - z * exp(-beta * energy)));

    return grandStatSum;
}

/** Calculates statistical sum in canonical ensemble based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> grandStatSum (grand statistical sum)
 * @param double phi (integration step)
 * @param double beta (1/T)
 * @param double saddlePoint
 * @return std::complex<double> - statistical sum
 * **/
inline complex<double> Droplet::statisticalSum(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                               const double &saddlePoint) const
{
    const complex<double> z = polar(saddlePoint, phi);
    return 1. / (2. * M_PI) * grandStatSum / pow(z, totalParticlesNumber + 1) * z;
}

/** This function calculates average stateOccupation of a specific state in canonical ensemble
 * based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> grandStatSum (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double foundedSaddlePoint (saddle point for the integration process)
 * @param int level (energy dropletEnergySpectrum level)
 * @return std::complex<double> - averageOccupation (average stateOccupation of a specific state)
 * **/
inline complex<double> Droplet::stateOccupation(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                                const double &saddlePoint, const unsigned &spectrumState) const
{
    const complex<double> z = polar(saddlePoint, phi);
    const double specificEnergyState = dropletEnergySpectrum[spectrumState];
    return 1. / (2. * M_PI) * z * exp(-beta * specificEnergyState) / (1. - z * exp(-beta * specificEnergyState)) *
           grandStatSum / pow(z, totalParticlesNumber + 1) * z;
}

/** This function calculates stateFluctuations of average stateOccupation of a specific state in canonical
 * ensemble based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> grandStatSum (grand partition function)
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature)
 * @param double foundedSaddlePoint (saddle point for the integration process)
 * @param int level (energy dropletEnergySpectrum level)
 * @return std::complex<double> - FL (stateFluctuations of a specific state)
 * **/
inline complex<double> Droplet::stateFluctuations(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                                  const double &saddlePoint, const unsigned &spectrumState) const
{
    const complex<double> z = polar(saddlePoint, phi);
    const double specificEnergyState = dropletEnergySpectrum[spectrumState];
    return 1. / (2. * M_PI) * z * exp(-beta * specificEnergyState) / (1. - z * exp(-beta * specificEnergyState)) *
           z * exp(-beta * specificEnergyState) / (1. - z * exp(-beta * specificEnergyState)) * grandStatSum / pow(z, totalParticlesNumber + 1) * z;
}

/** This function is used to calculate function to saddle point.
 * @param double z (searched saddle point)
 * @param double beta (beta coefficient related to temperature)
 * @param int Ntot (number of bosons)
 * @return double - function to saddle point
 * **/
inline double Droplet::saddlePointFunction(const double &z, const double &beta, const unsigned &Ntot) const
{
    double function = 0.; // Function on a saddle point

    for (const double &energy : dropletEnergySpectrum)
        function += z / (exp(beta * energy) - z);

    return function - (Ntot + 1);
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
 * @return double searchedSaddlePoint (saddle point for the integration process)
 * **/
double Droplet::saddlePoint(double start, double end, const double &beta, const unsigned &Ntot) const
{
    double searchedSaddlePoint = 0.; // Searched saddle point for the integration process

    if (saddlePointFunction(start, beta, Ntot) * saddlePointFunction(end, beta, Ntot) > 0.)
        return -1.;
    else
    {
        while ((end - start) * 0.5 > saddlePointAccuracy)
        {
            searchedSaddlePoint = (start + end) * 0.5;

            if (saddlePointFunction(searchedSaddlePoint, beta, Ntot) == 0)
                break;
            else if (saddlePointFunction(start, beta, Ntot) * saddlePointFunction(searchedSaddlePoint, beta, Ntot) < 0.)
                end = searchedSaddlePoint;
            else
                start = searchedSaddlePoint;
        }
    }

    return searchedSaddlePoint;
}

/** This function is used to calculate the specific state stateOccupation and
 * its basic properties. In this function we set beta coefficient and energy dropletEnergySpectrum level.
 * @param double beta (beta coefficient (beta = 1/T))
 * @param int bondState (specific state of energy dropletEnergySpectrum)
 * @return textfile (text file with properties)
 * **/
void Droplet::specificStateProperties(const double &T, const unsigned &spectrumState, const std::string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile << setprecision(5) << std::fixed << std::showpos;

    const double beta = 1. / T;
    const double foundedSaddlePoint = saddlePoint(0., 1., beta, totalParticlesNumber);
    
    double phi = 0.;
    complex<double> grandStatSum = 0.;

    double statSum = 0.;
    double averageOccupation = 0.;
    double fluctuations = 0.;

    #pragma omp parallel for private(phi, grandStatSum) reduction(+ : statSum, averageOccupation, fluctuations) num_threads(THREADS)
    for (unsigned k = 1; k <= integrationSteps; k++)
    {
        phi = integralLowerLimit + k * integrationAccuracy;
        grandStatSum = grandStatisticalSum(phi, beta, foundedSaddlePoint);
        statSum = statSum + statisticalSum(grandStatSum, phi, beta, foundedSaddlePoint).real();
        averageOccupation = averageOccupation + stateOccupation(grandStatSum, phi, beta, foundedSaddlePoint, spectrumState).real();
        fluctuations = fluctuations + stateFluctuations(grandStatSum, phi, beta, foundedSaddlePoint, spectrumState).real();
    }

    statSum = statSum * integrationAccuracy;
    averageOccupation = averageOccupation / statSum * integrationAccuracy;
    fluctuations = 2. / statSum * fluctuations * integrationAccuracy + averageOccupation - averageOccupation * averageOccupation;

    ofile << "Energy level:                                 " << spectrumState << "\n";
    ofile << "Total number of bosons:                       " << totalParticlesNumber << "\n";
    ofile << "Temperature:                                  " << T << "\n";
    ofile << "Saddle point:                                 " << foundedSaddlePoint << "\n";
    ofile << "Partition function:                           " << statSum << "\n";
    ofile << "Helmholtz free energy:                        " << -log(statSum) << "\n";
    ofile << "Average stateOccupation:                      " << averageOccupation << "\n";
    ofile << "Fluctuations of average stateOccupation:      " << fluctuations << "\n";
    ofile << "Dispersion of average stateOccupation:        " << sqrt(fluctuations) << "\n";
    ofile.close();
}

/** This function is used to check what happen with specific state
 * stateOccupation and stateFluctuations with temperature change.
 * @param int bondState (specific state of energy dropletEnergySpectrum)
 * @return Text file with Droplet properties
 * **/
void Droplet::specificStateTemperatureImpact(const unsigned &bondState, const double &Tp, const double &Tk, 
                                             const double &dT, const std::string &filename) const
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile << setprecision(5) << std::fixed << std::showpos;
    ofile << "T\t" << "F\t" << "<N>\t" << "Sigma\t\n"; 

    double phi = 0.;               
    complex<double> grandStatSum = 0.; 
    
    double statSum = 0.;           // Partition function
    double averageOccupation = 0.;           // Average stateOccupation
    double fluctuations = 0.;           // Fluctuations

    for (double T = Tp; T <= Tk; T += dT)
    {
        const double beta = 1. / T;                                                         
        const double foundedSaddlePoint = saddlePoint(0., 1., beta, totalParticlesNumber);

        #pragma omp parallel for private(phi, grandStatSum) reduction(+ : statSum, averageOccupation, fluctuations) num_threads(THREADS)
        for (unsigned k = 1; k <= integrationSteps; k++)
        {
            phi = integralLowerLimit + k * integrationAccuracy;
            grandStatSum = grandStatisticalSum(phi, beta, foundedSaddlePoint);
            statSum = statSum + statisticalSum(grandStatSum, phi, beta, foundedSaddlePoint).real();
            averageOccupation = averageOccupation + stateOccupation(grandStatSum, phi, beta, foundedSaddlePoint, bondState).real();
            fluctuations = fluctuations + stateFluctuations(grandStatSum, phi, beta, foundedSaddlePoint, bondState).real();
        }

        statSum = statSum * integrationAccuracy;
        averageOccupation = averageOccupation / statSum * integrationAccuracy;
        fluctuations = 2. / statSum * fluctuations * integrationAccuracy + averageOccupation - averageOccupation * averageOccupation;

        ofile << T << '\t' << -log(statSum) << '\t' << averageOccupation << '\t' << fluctuations << '\n';
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
std::stringstream Droplet::dropletWidth(const double &T, double &particlesInDroplet)
{
    // Counting bond levels
    unsigned maxBondStates = 0;
    while (dropletEnergySpectrum[maxBondStates] - systemPotential < saddlePointAccuracy)
        maxBondStates++;

    const double beta = 1. / T;
    const double foundedSaddlePoint = saddlePoint(0., 1., beta, totalParticlesNumber);
    double helmholtzFreeEnergy = 0.;
    particlesInDroplet = 0.;

    double phi = 0.;
    complex<double> grandStatSum = 0.;

    double statSum = 0.;           
    double averageOccupation = 0.;                    

    for (unsigned bondState = 0; bondState < maxBondStates; bondState++)
    {
        #pragma omp parallel for private(phi, grandStatSum) reduction(+ : statSum, averageOccupation) num_threads(THREADS)
        for (unsigned k = 1; k <= integrationSteps; k++)
        {
            phi = integralLowerLimit + k * integrationAccuracy;
            grandStatSum = grandStatisticalSum(phi, beta, foundedSaddlePoint);
            statSum = statSum + statisticalSum(grandStatSum, phi, beta, foundedSaddlePoint).real();
            averageOccupation = averageOccupation + stateOccupation(grandStatSum, phi, beta, foundedSaddlePoint, bondState).real();
        }

        statSum = statSum * integrationAccuracy;
        averageOccupation = averageOccupation / statSum * integrationAccuracy;

        particlesInDroplet += averageOccupation;
        helmholtzFreeEnergy += -log(statSum);
    }

    stringstream output;
    output << T << "\t" << scaleCoefficient * particlesInDroplet << "\t" << helmholtzFreeEnergy << "\n";
    return output;
}