#ifndef DROPLET_H
#define DROPLET_H

#include <vector>
#include <complex>
#include <sstream>

class Droplet
{
private:
    unsigned integrationSteps;                 // Integration steps
    double integralLowerLimit;                 // Lower limit of the integral
    double integralUpperLimit;                 // Upper limit of the integral
    double integrationAccuracy;                // Integration step
    unsigned Ntot;                             // Total number of bosons
    double systemPotential;                    // Potential
    double width;                              // Initial Droplet width
    double scaleCoefficient;                   // Constant scale coefficient beetwen Ntot and width
    std::vector<double> dropletEnergySpectrum; // Energy spectrum, first state is a Droplet state
    double saddlePointAccuracy;                // Saddle point accuracy

    std::complex<double> grandStatisticalSum(const double &phi, const double &beta, const double &saddlePoint) const;

    std::complex<double> statisticalSum(const std::complex<double> &gSS, const double &phi, const double &beta,
                                        const double &saddlePoint) const;

    std::complex<double> stateOccupation(const std::complex<double> &gSS, const double &phi, const double &beta,
                                         const double &saddlePoint, const unsigned &specState) const;

    std::complex<double> stateFluctuations(const std::complex<double> &gSS, const double &phi, const double &beta,
                                           const double &saddlePoint, const unsigned &specState) const;

    double saddlePointFunction(const double &z, const double &beta, const unsigned &Ntot) const;
    double saddlePoint(double start, double end, const double &beta, const unsigned &Ntot) const;

public:
    Droplet(const unsigned &NTot, const double &Width, const double &V);
    void setSpectrum(const std::vector<double> &Spectrum);
    void stateProperties(const double &T, const unsigned &specState, std::string filename);
    void stateTemperatureImpact(const unsigned &specState, const double &Tp, const double &Tk, const double &dT, std::string filename);
    std::stringstream dropletWidth(const double &T, double &NDroplet);
};

#endif // DROPLET_H