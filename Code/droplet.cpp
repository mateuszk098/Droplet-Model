#define _USE_MATH_DEFINES
#define THREADS 12

#include "Droplet.h"
#include <omp.h>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

Droplet::Droplet(const unsigned &totalParticlesNumber, const double &dropletWidth, const double &systemPotential) : 
totalParticlesNumber(totalParticlesNumber), dropletWidth(dropletWidth), systemPotential(systemPotential),
integralLowerLimit(-M_PI), integralUpperLimit(M_PI), saddlePointAccuracy(1e-6)
{
    scaleCoefficient = dropletWidth / totalParticlesNumber;
    integrationSteps = totalParticlesNumber * 100;
    integrationAccuracy = (integralUpperLimit - integralLowerLimit) / integrationSteps;
}

void Droplet::setSpectrum(const std::vector<double> &dropletEnergySpectrum) { this->dropletEnergySpectrum = dropletEnergySpectrum; }

/** 
 * Calculates grand statistical sum in grand canonical ensemble
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature beta = 1/T)
 * @param double saddlePoint
 * @return std::complex<double> - grand statistical sum
 */
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

/** 
 * Calculates statistical sum in canonical ensemble based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> grand statistical sum
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature beta = 1/T)
 * @param double saddlePoint
 * @return std::complex<double> - statistical sum
 */
inline complex<double> Droplet::statisticalSum(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                               const double &saddlePoint) const
{
    const complex<double> z = polar(saddlePoint, phi);
    return 1. / (2. * M_PI) * grandStatSum / pow(z, totalParticlesNumber + 1) * z;
}

/** 
 * Calculates average occupation of a specific state in canonical ensemble
 * based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> grand statistical sum
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature beta = 1/T)
 * @param double saddlePoint
 * @param unsigned spectrum specific state
 * @return std::complex<double> - average occupation of a specific state
 */
inline complex<double> Droplet::stateOccupation(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                                const double &saddlePoint, const unsigned &spectrumState) const
{
    const complex<double> z = polar(saddlePoint, phi);
    const double specificEnergyState = dropletEnergySpectrum[spectrumState];
    return 1. / (2. * M_PI) * z * exp(-beta * specificEnergyState) / (1. - z * exp(-beta * specificEnergyState)) *
           grandStatSum / pow(z, totalParticlesNumber + 1) * z;
}

/** 
 * Calculates fluctuations of a specific state in canonical
 * ensemble based on grand statistical sum using path integral forumaltion.
 * @param std::complex<double> grand statistical sum
 * @param double phi (integration step)
 * @param double beta (beta coefficient related to temperature beta = 1/T)
 * @param double saddlePoint
 * @param unsigned spectrum specific state
 * @return std::complex<double> - fluctuations of a specific state
 */
inline complex<double> Droplet::stateFluctuations(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                                  const double &saddlePoint, const unsigned &spectrumState) const
{
    const complex<double> z = polar(saddlePoint, phi);
    const double specificEnergyState = dropletEnergySpectrum[spectrumState];
    return 1. / (2. * M_PI) * z * exp(-beta * specificEnergyState) / (1. - z * exp(-beta * specificEnergyState)) *
           z * exp(-beta * specificEnergyState) / (1. - z * exp(-beta * specificEnergyState)) * 
           grandStatSum / pow(z, totalParticlesNumber + 1) * z;
}

/**
 * Calculates function to saddle point.
 * @param double z (searched saddle point)
 * @param double beta (beta coefficient related to temperature beta = 1/T)
 * @param int totalParticlesNumber
 * @return double - function to saddle point
 */
inline double Droplet::saddlePointFunction(const double &z, const double &beta, const unsigned &totalParticlesNumber) const
{
    double function = 0.;

    for (const double &energy : dropletEnergySpectrum)
        function += z / (exp(beta * energy) - z);

    return function - (totalParticlesNumber + 1);
}

/** 
 * Bisection method to find saddle point to integration process.
 * @param double start (start point of a range)
 * @param double end (end point of a range)
 * @param double beta (beta coefficient related to temperature beta = 1/T)
 * @param int totalParticlesNumber
 * @return double - searchedSaddlePoint
 */
double Droplet::saddlePoint(double start, double end, const double &beta, const unsigned &totalParticlesNumber) const
{
    double searchedSaddlePoint = 0.;

    if (saddlePointFunction(start, beta, totalParticlesNumber) * saddlePointFunction(end, beta, totalParticlesNumber) > 0.)
        return -1.;
    else
    {
        while ((end - start) * 0.5 > saddlePointAccuracy)
        {
            searchedSaddlePoint = (start + end) * 0.5;

            if (saddlePointFunction(searchedSaddlePoint, beta, totalParticlesNumber) == 0)
                break;
            else if (saddlePointFunction(start, beta, totalParticlesNumber) * saddlePointFunction(searchedSaddlePoint, beta, totalParticlesNumber) < 0.)
                end = searchedSaddlePoint;
            else
                start = searchedSaddlePoint;
        }
    }

    return searchedSaddlePoint;
}

/** 
 * Calculates the specific state occupation and its basic properties. 
 * @param double temperature
 * @param unsigned spectrum specific state
 * @param string filename where to save results
 * @return Text file with droplet properties
 */
void Droplet::specificStateProperties(const double &T, const unsigned &spectrumState, const std::string &filename)
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile << setprecision(5) << std::fixed << std::showpos;

    const double beta = 1. / T;
    const double foundSaddlePoint = saddlePoint(0., 1., beta, totalParticlesNumber);

    double phi = 0.;
    complex<double> grandStatSum = 0.;

    double statSum = 0.;
    double averageOccupation = 0.;
    double fluctuations = 0.;

    #pragma omp parallel for private(phi, grandStatSum) reduction(+ : statSum, averageOccupation, fluctuations) num_threads(THREADS)
    for (unsigned k = 1; k <= integrationSteps; k++)
    {
        phi = integralLowerLimit + k * integrationAccuracy;
        grandStatSum = grandStatisticalSum(phi, beta, foundSaddlePoint);
        statSum = statSum + statisticalSum(grandStatSum, phi, beta, foundSaddlePoint).real();
        averageOccupation = averageOccupation + stateOccupation(grandStatSum, phi, beta, foundSaddlePoint, spectrumState).real();
        fluctuations = fluctuations + stateFluctuations(grandStatSum, phi, beta, foundSaddlePoint, spectrumState).real();
    }

    statSum = statSum * integrationAccuracy;
    averageOccupation = averageOccupation / statSum * integrationAccuracy;
    fluctuations = 2. / statSum * fluctuations * integrationAccuracy + averageOccupation - averageOccupation * averageOccupation;

    ofile << "Energy level:                                 " << spectrumState << "\n";
    ofile << "Total number of bosons:                       " << totalParticlesNumber << "\n";
    ofile << "Temperature:                                  " << T << "\n";
    ofile << "Saddle point:                                 " << foundSaddlePoint << "\n";
    ofile << "Partition function:                           " << statSum << "\n";
    ofile << "Helmholtz free energy:                        " << -log(statSum) << "\n";
    ofile << "Average stateOccupation:                      " << averageOccupation << "\n";
    ofile << "Fluctuations of average stateOccupation:      " << fluctuations << "\n";
    ofile << "Dispersion of average stateOccupation:        " << sqrt(fluctuations) << "\n";
    ofile.close();
}

/** 
 * Calculates what happen with specific state with temperature change.
 * @param unsigned specific spectrum state
 * @param double start temperature
 * @param double end temperature
 * @param double temperature step
 * @param string filename where to save results
 * @return Text file with Droplet properties
 */
void Droplet::specificStateTemperatureImpact(const unsigned &spectrumState, const double &startTemperature, const double &endTemperature,
                                             const double &temperatureStep, const std::string &filename) const
{
    ofstream ofile(filename.c_str(), ios::out);
    ofile << setprecision(5) << std::fixed << std::showpos;
    ofile << "T\t" << "F\t" << "<N>\t" << "Sigma\t\n";

    double phi = 0.;
    complex<double> grandStatSum = 0.;

    double statSum = 0.;           
    double averageOccupation = 0.; 
    double fluctuations = 0.;      

    for (double T = startTemperature; T <= endTemperature; T += temperatureStep)
    {
        const double beta = 1. / T;
        const double foundSaddlePoint = saddlePoint(0., 1., beta, totalParticlesNumber);

        #pragma omp parallel for private(phi, grandStatSum) reduction(+ : statSum, averageOccupation, fluctuations) num_threads(THREADS)
        for (unsigned k = 1; k <= integrationSteps; k++)
        {
            phi = integralLowerLimit + k * integrationAccuracy;
            grandStatSum = grandStatisticalSum(phi, beta, foundSaddlePoint);
            statSum = statSum + statisticalSum(grandStatSum, phi, beta, foundSaddlePoint).real();
            averageOccupation = averageOccupation + stateOccupation(grandStatSum, phi, beta, foundSaddlePoint, spectrumState).real();
            fluctuations = fluctuations + stateFluctuations(grandStatSum, phi, beta, foundSaddlePoint, spectrumState).real();
        }

        statSum = statSum * integrationAccuracy;
        averageOccupation = averageOccupation / statSum * integrationAccuracy;
        fluctuations = 2. / statSum * fluctuations * integrationAccuracy + averageOccupation - averageOccupation * averageOccupation;

        ofile << T << '\t' << -log(statSum) << '\t' << averageOccupation << '\t' << fluctuations << '\n';
    }

    ofile.close();
}

/** 
 * Calculate droplet width in the specific temperature.
 * Quantum Droplet consist of first droplet state and every bond states.
 * @param double temperature
 * @param double particlesInDroplet
 * @return Text file with droplet properties
 */
std::tuple<double, double> Droplet::calculateDropletWidth(const double &T, double &particlesInDroplet)
{
    // Counting bond levels
    unsigned maxBondStates = 0;
    while (dropletEnergySpectrum[maxBondStates] - systemPotential < 1e-6)
        maxBondStates++;

    const double beta = 1. / T;
    const double foundSaddlePoint = saddlePoint(0., 1., beta, totalParticlesNumber);
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
            grandStatSum = grandStatisticalSum(phi, beta, foundSaddlePoint);
            statSum = statSum + statisticalSum(grandStatSum, phi, beta, foundSaddlePoint).real();
            averageOccupation = averageOccupation + stateOccupation(grandStatSum, phi, beta, foundSaddlePoint, bondState).real();
        }

        statSum = statSum * integrationAccuracy;
        averageOccupation = averageOccupation / statSum * integrationAccuracy;

        particlesInDroplet += averageOccupation;
        helmholtzFreeEnergy += -log(statSum);
    }

    return make_tuple(scaleCoefficient * particlesInDroplet, helmholtzFreeEnergy);
}