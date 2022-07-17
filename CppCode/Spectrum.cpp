#define _USE_MATH_DEFINES

#include "Spectrum.h"
#include <algorithm>
#include <cmath>
using namespace std;

Spectrum::Spectrum(const double &bisectionAccuracy, const unsigned &subBisectionNumber) : bisectionAccuracy(bisectionAccuracy),
                                                                                          subBisectionNumber(subBisectionNumber) {}

inline void Spectrum::setSubBisectionNumber(const unsigned &subBisectionNumber) { this->subBisectionNumber = subBisectionNumber; }

/** Returns function to roots of even bond solutions
 * @param double x (function variable)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of even bound solutions)
 * **/
inline double Spectrum::functionToEvenBondRoots(const double &x, const double &L, const double &a, const double &c) const
{
    return x / tanh((1. - (L / a)) * x) + sqrt(c - x * x) * tan(sqrt(c - x * x));
}

/** Returns function to roots of odd bound solutions
 * @param double x (function variable)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of odd bound solutions)
 * **/
inline double Spectrum::functionToOddBondRoots(const double &x, const double &L, const double &a, const double &c) const
{
    return x / tanh((1. - (L / a)) * x) - sqrt(c - x * x) / tan(sqrt(c - x * x));
}

/** This method returns function to roots of even scatter solutions
 * @param double x (function variable)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of even scatter solutions)
 * **/
inline double Spectrum::functionToEvenScatterRoots(const double &x, const double &L, const double &a, const double &c) const
{
    return x / tan((1. - (L / a)) * x) + sqrt(c + x * x) * tan(sqrt(c + x * x));
}

/** This method returns function to roots of odd scatter solutions
 * @param double x (function variable)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of odd scatter solutions)
 * **/
inline double Spectrum::functionToOddScatterRoots(const double &x, const double &L, const double &a, const double &c) const
{
    return x / tan((1 - (L / a)) * x) - sqrt(c + x * x) / tan(sqrt(c + x * x));
}

/** Uses bisection method to find root of function in specific subrange (xp; xk)
 * @param double xp (begin of the range)
 * @param double xk (end of the range)
 * @param function fptr (pointer to appropriate function)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double root (root of function in the range (xp; xk), if method not found root in the range, then -1 is returned)
 * **/
double Spectrum::findFunctionRoot(double xp, double xk, function fptr, const double &L, const double &a, const double &c) const
{
    double functionRoot = 0.;

    // We are looking for the root in range (xp, xk)
    // Bisection method assumes root not exist if f(tp) * f(tk) > 0
    if ((this->*fptr)(xp, L, a, c) * (this->*fptr)(xk, L, a, c) > 0.)
        return -1;
    // Else possibly we have root there
    else
    {
        while ((xk - xp) * 0.5 > bisectionAccuracy)
        {
            functionRoot = (xp + xk) * 0.5;

            if ((this->*fptr)(functionRoot, L, a, c) == 0)
                break;
            else if ((this->*fptr)(xp, L, a, c) * (this->*fptr)(functionRoot, L, a, c) < 0)
                xk = functionRoot;
            else
                xp = functionRoot;
        }
    }

    return functionRoot;
}

/** This method using multi bisection method to find every roots of function in whole range (xp; xk)
 * @param double xp (begin of the range)
 * @param double xk (end of the range)
 * @param function fptr (pointer to appropriate function)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return std::vector<double> allFunctionRoots (roots of function in the whole range (xp; xk))
 * **/
vector<double> Spectrum::findAllFunctionRoots(double xp, double xk, function fptr, const double &L, const double &a, const double &c) const
{
    const double rangeInterval = abs(xk - xp) / subBisectionNumber;
    vector<double> allFunctionRoots;
    double functionRoot = 0.;

    for (unsigned i = 0; i < subBisectionNumber; i++)
    {
        functionRoot = findFunctionRoot(xp, xp + rangeInterval, fptr, L, a, c);
        // Omit -1 roots might be weird but we know we are looking for only positive roots so this is a correct approach
        // Apart from that, function find_root must returns something. So I adopt if root does not exist, -1 will be returned
        // But the rangeInterval must be increased
        if (functionRoot == -1)
        {
            xp += rangeInterval;
            continue;
        }
        allFunctionRoots.push_back(functionRoot);
        xp += rangeInterval;
    }

    return allFunctionRoots;
}

/** Find solutions of finite well potential problem
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double V (depth of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @param double xp (begin of the range)
 * @param double xk (end of the range)
 * @return std::vector<double> solutions of finite well potential problem
 * **/
vector<double> Spectrum::calculateWellSpectrum(const double &L, const double &a, const double &V, const double &c,
                                               const double &xp, const double &xkb, const double &xks) const
{
    vector<double> evenBondRoots(findAllFunctionRoots(xp, xkb, functionToEvenBondRoots, L, a, c));
    vector<double> oddBondRoots(findAllFunctionRoots(xp, xkb, functionToOddBondRoots, L, a, c));
    vector<double> evenScatterRoots(findAllFunctionRoots(xp, xks, functionToEvenScatterRoots, L, a, c));
    vector<double> oddScatterRoots(findAllFunctionRoots(xp, xks, functionToOddScatterRoots, L, a, c));

    // I have to remove roots that have been incorrectly identified as roots, because they are asymptotes
    // I just check if function in this root is greater than 1, if so, it is rather the asymptote
    // I will use appropriate lambda functions
    const auto lambdaEbr = [this, L, a, c](const double &x){ return abs(functionToEvenBondRoots(x, L, a, c)) > 1.0; };
    const auto lambdaObr = [this, L, a, c](const double &x){ return abs(functionToOddBondRoots(x, L, a, c)) > 1.0; };
    const auto lambdaEsr = [this, L, a, c](const double &x){ return abs(functionToEvenScatterRoots(x, L, a, c)) > 1.0; };
    const auto lambdaOsr = [this, L, a, c](const double &x){ return abs(functionToOddScatterRoots(x, L, a, c)) > 1.0; };

    // Now I remove these asymptotes from appropriate vectors
    evenBondRoots.erase(remove_if(evenBondRoots.begin(), evenBondRoots.end(), lambdaEbr), evenBondRoots.end());
    oddBondRoots.erase(remove_if(oddBondRoots.begin(), oddBondRoots.end(), lambdaObr), oddBondRoots.end());
    evenScatterRoots.erase(remove_if(evenScatterRoots.begin(), evenScatterRoots.end(), lambdaEsr), evenScatterRoots.end());
    oddScatterRoots.erase(remove_if(oddScatterRoots.begin(), oddScatterRoots.end(), lambdaOsr), oddScatterRoots.end());

    // I combine bond solutions and scatter solutions
    vector<double> bondRoots;
    vector<double> scatterRoots;
    bondRoots.insert(bondRoots.begin(), evenBondRoots.begin(), evenBondRoots.end());
    bondRoots.insert(bondRoots.end(), oddBondRoots.begin(), oddBondRoots.end());
    scatterRoots.insert(scatterRoots.begin(), evenScatterRoots.begin(), evenScatterRoots.end());
    scatterRoots.insert(scatterRoots.end(), oddScatterRoots.begin(), oddScatterRoots.end());

    // I evaluate bond energy levels
    const auto lambdaBondEnergy = [a](const double &x){ return -(x * x) / (2. * a * a); };
    transform(bondRoots.begin(), bondRoots.end(), bondRoots.begin(), lambdaBondEnergy);

    // I evaluate scatter energy levels
    const auto lambdaScatterEnergy = [a](const double &x){ return (x * x) / (2. * a * a); };
    transform(scatterRoots.begin(), scatterRoots.end(), scatterRoots.begin(), lambdaScatterEnergy);

    // I combine bound and scatter spectrum and sort
    vector<double> full_spectrum;
    full_spectrum.insert(full_spectrum.begin(), bondRoots.begin(), bondRoots.end());
    full_spectrum.insert(full_spectrum.end(), scatterRoots.begin(), scatterRoots.end());
    sort(full_spectrum.begin(), full_spectrum.end());

    // Eventually I translate a whole spectrum by potential depth
    const auto lambdaTranslateSpectrum = [V](const double &x){ return x + V; };
    transform(full_spectrum.begin(), full_spectrum.end(), full_spectrum.begin(), lambdaTranslateSpectrum);

    // And I add first droplet energy state equal zero
    full_spectrum.insert(full_spectrum.begin(), 0.);

    return full_spectrum;
}

/** Returns spectrum from master thesis
 * @param double N (number of partciles whose create droplet)
 * @param double gdd (dipolar coupling constant)
 * @param double V (chemical potential -> potential)
 * @param double L (width of the system with surroundings)
 * @param unsigned scatterStatesNumber 
 * @return std::vector<double> spectrum from master thesis
 * **/
std::vector<double> Spectrum::calculateMasterThesisSpectrum(const double &N, const double &gdd, const double &V, const double &L,
                                                            const unsigned &scatterStatesNumber) const
{
    unsigned m = 1;
    double E = 0.;
    std::vector<double> masterThesisSpectrum;

    // First bound level is a droplet level
    masterThesisSpectrum.push_back(0.);

    // Bound levels
    while (true)
    {
        E = (m / N) * (9. / (4. * M_PI * M_PI)) * sqrt(0.333333 + 0.25 * (m * m) / (N * N)) * gdd * gdd;

        if (E - V < 1e-6)
        {
            masterThesisSpectrum.push_back(E);
            m++;
        }
        else
        {    
            break;
        }
    }

    // Scatter levels
    for (unsigned k = 1; k <= scatterStatesNumber; k++)
    {
        // PBC
        E = 0.5 * (2. * M_PI * k / L) * (2. * M_PI * k / L) + V;
        masterThesisSpectrum.push_back(E);
        masterThesisSpectrum.push_back(E);

        // OBC
        // E = 0.5 * (M_PI * k / L) * (M_PI * k / L) + V;
        // full_spectrum.push_back(E);
    }

    return masterThesisSpectrum;
}
