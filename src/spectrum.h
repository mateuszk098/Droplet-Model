#ifndef SPECTRUM
#define SPECTRUM
#include <vector>

// This class is used to calculate energy spectrum in finite well potential problem
class Spectrum
{
private:
    // Accuracy of the bisection method
    double bisectionAccuracy;
    unsigned subBisectionNumber;

    typedef double (Spectrum::*function)(const double &x, const double &L, const double &a, const double &cons) const;

    double functionToEvenBondRoots(const double &x, const double &L, const double &a, const double &c) const;
    double functionToOddBondRoots(const double &x, const double &L, const double &a, const double &c) const;
    double functionToEvenScatterRoots(const double &x, const double &L, const double &a, const double &c) const;
    double functionToOddScatterRoots(const double &x, const double &L, const double &a, const double &c) const;

    double findFunctionRoot(double xp, double xk, function fptr, const double &L, const double &a, const double &c) const;
    std::vector<double> findAllFunctionRoots(double xp, double end, function fptr, const double &L, const double &a, const double &c) const;

public:
    Spectrum(const double &bisectionAccuracy, const unsigned &subBisectionNumber);
    void setSubBisectionNumber(const unsigned &subBisectionNumber);
    std::vector<double> calculateWellSpectrum(const double &L, const double &a, const double &V, const double &c,
                                              const double &xp, const double &xkb, const double &xks) const;
    std::vector<double> calculateMasterThesisSpectrum(const double &N, const double &gdd, const double &V, const double &L,
                                                      const unsigned &scatterStatesNumber) const;
};

#endif // WELL_SPECTRUM_H