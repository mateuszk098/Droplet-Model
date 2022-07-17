#ifndef DROPLET_H
#define DROPLET_H

#include <vector>
#include <complex>
#include <sstream>

class Droplet
{
private:
    unsigned totalParticlesNumber;
    double dropletWidth;
    double systemPotential;

    double integralLowerLimit;
    double integralUpperLimit;
    double saddlePointAccuracy;
    double scaleCoefficient;
    unsigned integrationSteps;
    double integrationAccuracy;

    std::vector<double> dropletEnergySpectrum;

    std::complex<double> grandStatisticalSum(const double &phi, const double &beta, const double &saddlePoint) const;

    std::complex<double> statisticalSum(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                        const double &saddlePoint) const;

    std::complex<double> stateOccupation(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                         const double &saddlePoint, const unsigned &spectrumState) const;

    std::complex<double> stateFluctuations(const std::complex<double> &grandStatSum, const double &phi, const double &beta,
                                           const double &saddlePoint, const unsigned &spectrumState) const;

    double saddlePointFunction(const double &z, const double &beta, const unsigned &Ntot) const;
    double saddlePoint(double start, double end, const double &beta, const unsigned &Ntot) const;

public:
    Droplet(const unsigned &totalParticlesNumber, const double &dropletWidth, const double &systemPotential);

    void setSpectrum(const std::vector<double> &dropletEnergySpectrum);

    void specificStateProperties(const double &T, const unsigned &spectrumState, const std::string &filename);

    void specificStateTemperatureImpact(const unsigned &specState, const double &Tp, const double &Tk,
                                        const double &dT, const std::string &filename) const;

    std::stringstream calculateDropletWidth(const double &T, double &particlesInDroplet);
};

#endif // DROPLET_H