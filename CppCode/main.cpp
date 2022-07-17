/** Author: Mateusz Kowalczyk
 * This program calculates interesting properties of a new quantum state,
 * quantum droplet. The main method of calculating is based on performing
 * the calculations of canonical ensemble based on grand canonical ensemble.
 * For this purpose, I use the outputPath integration formula based on Cauchy's
 * integration formula. This program using multithreads.
 **/

#define _USE_MATH_DEFINES

#include "Droplet.h"
#include "Spectrum.h"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

void loadParameters(string &inputFilename, string &outputFilename, unsigned &spectrumType, unsigned &totalParticlesNumber,
                    unsigned &spectrumEnergyState, double &gdd, double &startTemperature, double &endTemperature,
                    double &temperatureStep, unsigned &thermalisationSteps)
{
    string line, tmp;
    ifstream input(inputFilename, ios::in);
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> tmp >> outputFilename >> tmp >> spectrumType >> tmp >> totalParticlesNumber >> tmp >> spectrumEnergyState;
            input >> tmp >> gdd >> tmp >> startTemperature >> tmp >> endTemperature >> tmp >> temperatureStep >> tmp >> thermalisationSteps;
        }
        input.close();
    }
    else
        cerr << "Cannot open file!" << '\n';
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Too few arguments!" << '\n';
        return EXIT_FAILURE;
    }
    ios::sync_with_stdio(0);
    cin.tie(0);

    string inputFilename = argv[1];
    string outputFilename = "";

    // Droplet parameters
    unsigned spectrumType = 0;
    unsigned totalParticlesNumber = 0;
    unsigned spectrumEnergyState = 0;
    double gdd = 0.;
    double startTemperature = 0.;
    double endTemperature = 0.;
    double temperatureStep = 0.;
    unsigned thermalisationSteps = 0;

    loadParameters(inputFilename, outputFilename, spectrumType, totalParticlesNumber, spectrumEnergyState,
                   gdd, startTemperature, endTemperature, temperatureStep, thermalisationSteps);

    double particlesInDroplet = static_cast<double>(totalParticlesNumber);

    // Model parameters
    double a = 2. * M_PI * M_PI * totalParticlesNumber / (3. * gdd); // Initial droplet width
    double V = 3. / (8. * M_PI * M_PI) * gdd * gdd;                  // Potential
    double L = 10. * a;                                              // Walls
    double c = 2. * V * a * a;                                       // Constant coefficient in model

    // Parameters used in finite potential well
    const unsigned subBisectionNumber = 20000;
    const double bisectionAccuracy = 1e-6;
    const double xp = bisectionAccuracy;             // Begin of range
    const double xks = 100. - bisectionAccuracy;     // End of range for scatter spectrumEnergyStates
    double xkb = sqrt(c * 0.25) - bisectionAccuracy; // End of range for bond spectrumEnergyStates

    // Initialise objects
    vector<double> full_spectrum;
    Droplet dropletInstance(totalParticlesNumber, a, V);
    Spectrum spectrumInstance(bisectionAccuracy, subBisectionNumber);

    // ----------------------- SIMULATION LOOP ----------------------------------------------
    const string outputPath = "../Data/" + outputFilename;
    ofstream outputFile(outputPath.c_str(), ios::out);

    outputFile << setprecision(5);
    auto tp = high_resolution_clock::now();

    // Main temperature loop
    for (double T = startTemperature; T < endTemperature; T += temperatureStep)
    {
        for (unsigned index = 0; index < thermalisationSteps; index++)
        {
            // -------------- POTENTIAL WELL SPECTRUM ---------------------------------------
            // Parameters a, c and xkb change with every step
            if (spectrumType == 1)
            {
                a = 2.0 * M_PI * M_PI * particlesInDroplet / (3.0 * gdd);
                c = 2.0 * V * a * a;
                xkb = sqrt(c * 0.25) - bisectionAccuracy;
                full_spectrum = spectrumInstance.calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks);
            }
            // -------------- END OF POTENTIAL WELL SPECTRUM --------------------------------

            // -------------- REAL SPECTRUM FROM MASTER THESIS ------------------------------
            else if (spectrumType == 2)
            {
                unsigned spectrumEnergyStates = 1000;
                full_spectrum = spectrumInstance.calculateMasterThesisSpectrum(particlesInDroplet, gdd, V, L, spectrumEnergyStates);
            }
            // -------------- END OF REAL SPECTRUM FROM MASTER THESIS -----------------------

            dropletInstance.setSpectrum(full_spectrum);
            double dropletWidth, helmholtzFreeEnergy;
            tie(dropletWidth, helmholtzFreeEnergy) = dropletInstance.calculateDropletWidth(T, particlesInDroplet);

            if (index == thermalisationSteps - 1)
                outputFile << T << '\t' << dropletWidth << '\t' << helmholtzFreeEnergy << '\n';
        }
    }

    outputFile.close();
    auto tk = high_resolution_clock::now();           // Execution time - stop
    duration<double, std::milli> ms_double = tk - tp; // Getting  milliseconds as a double
    cout << "Execution time on the CPU: " << ms_double.count() * 1e-3 << " s";
    // ----------------------- END OF SIMULATION LOOP ---------------------------------------

    // ----------------------- TEST OTHER FUNCTIONS -----------------------------------------
    // xkb = sqrt(c * 0.25) - bisectionAccuracy;
    // full_spectrum = spectrumInstance.calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks); // Quantum well
    // full_spectrum = rs->spectrum(particlesInDroplet, gdd, V, L, 1000); // Real
    // for (int i = 0; i < 1000; i++) // 1D condensate spectrum
    //     full_spectrum.push_back(i);

    // dropletInstance.setSpectrum(full_spectrum);
    // dropletInstance.specificStateProperties(Tp, spectrumEnergyState, outputFilename);
    // dropletInstance.specificStateTemperatureImpact(spectrumEnergyState, Tp, Tk, dT, outputFilename);
    // ----------------------- END OF TEST OTHER FUNCTIONS ----------------------------------

    return EXIT_SUCCESS;
}
