/**
 * @author: Mateusz Kowalczyk
 * @name: Droplet Model
 * This program calculates interesting properties of a new quantum state,
 * quantum droplet. The main method of calculating is based on performing
 * the calculations of canonical ensemble based on grand canonical ensemble.
 * For this purpose, I use the path integration formula based on Cauchy's
 * integration formula. This program using multithreads.
 * Compilation: c++ @flags
 * Usage: ./main parameters.txt
 */

#define _USE_MATH_DEFINES

#include "droplet.h"
#include "spectrum.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

/*
Description of some parameters.
- option:
    1 - Investigation of the whole quantum droplet model (temperature influence to droplet width).
    2 - Investigation of some interesting properties of exactly one quantum state of the droplet.
    3 - Investigation of the temperature influence to the number of particles in exactly one state.
- spectrumType:
    1 - Quantum Well
    2 - Spectrum from Master Thesis
- thermalisationSteps - after how many steps in for `temperatureStep` save the result
    3 - or more, the safest tested is 5.
*/
void loadParameters(string &inputFilename, string &outputFilename, unsigned &option, unsigned &spectrumType,
                    unsigned &thermalisationSteps, unsigned &totalParticlesNumber, unsigned &spectrumEnergyState,
                    double &gdd, double &startTemperature, double &endTemperature, double &temperatureStep)
{
    string line, tmp;
    ifstream input(inputFilename, ios::in);
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> tmp >> outputFilename >> tmp >> option >> tmp >> spectrumType >> tmp >> thermalisationSteps;
            input >> tmp >> totalParticlesNumber >> tmp >> spectrumEnergyState >> tmp >> gdd;
            input >> tmp >> startTemperature >> tmp >> endTemperature >> tmp >> temperatureStep >> tmp;
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
    unsigned option = 0;
    unsigned spectrumType = 0;
    unsigned thermalisationSteps = 0;
    unsigned totalParticlesNumber = 0;
    unsigned spectrumEnergyState = 0;
    double gdd = 0.;
    double startTemperature = 0.;
    double endTemperature = 0.;
    double temperatureStep = 0.;

    loadParameters(inputFilename, outputFilename, option, spectrumType, thermalisationSteps, totalParticlesNumber,
                   spectrumEnergyState, gdd, startTemperature, endTemperature, temperatureStep);

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
    const double xks = 120.;                         // End of range for scatter spectrumEnergyStates
    double xkb = sqrt(c * 0.25) - bisectionAccuracy; // End of range for bond spectrumEnergyStates

    Droplet dropletInstance(totalParticlesNumber, a, V);
    Spectrum spectrumInstance(bisectionAccuracy, subBisectionNumber);

    // ----------------------- SIMULATION LOOP ----------------------------------------------
    const string outputPath = "../data/" + outputFilename;
    ofstream outputFile(outputPath.c_str(), ios::out);
    outputFile << setprecision(5);

    auto tp = high_resolution_clock::now();

    switch (option)
    {
    case 1:
    {
        switch (spectrumType)
        {
        case 1:
        {
            for (double T = startTemperature; T < endTemperature; T += temperatureStep)
            {
                for (unsigned index = 0; index < thermalisationSteps; index++)
                {
                    a = 2.0 * M_PI * M_PI * particlesInDroplet / (3.0 * gdd);
                    c = 2.0 * V * a * a;
                    xkb = sqrt(c * 0.25) - bisectionAccuracy;
                    const vector<double> dropletSpectrum(spectrumInstance.calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks));

                    dropletInstance.setSpectrum(dropletSpectrum);
                    double dropletWidth, helmholtzFreeEnergy;
                    tie(dropletWidth, helmholtzFreeEnergy) = dropletInstance.calculateDropletWidth(T, particlesInDroplet);

                    if (index == thermalisationSteps - 1)
                        outputFile << T << '\t' << dropletWidth << '\t' << helmholtzFreeEnergy << '\n';
                }
            }
            break;
        }
        case 2:
        {
            for (double T = startTemperature; T < endTemperature; T += temperatureStep)
            {
                for (unsigned index = 0; index < thermalisationSteps; index++)
                {
                    unsigned spectrumEnergyStates = 500;
                    const vector<double> dropletSpectrum(spectrumInstance.calculateMasterThesisSpectrum(particlesInDroplet, gdd, V, L, spectrumEnergyStates));

                    dropletInstance.setSpectrum(dropletSpectrum);
                    double dropletWidth, helmholtzFreeEnergy;
                    tie(dropletWidth, helmholtzFreeEnergy) = dropletInstance.calculateDropletWidth(T, particlesInDroplet);

                    if (index == thermalisationSteps - 1)
                        outputFile << T << '\t' << dropletWidth << '\t' << helmholtzFreeEnergy << '\n';
                }
            }
            break;
        }

        default:
        {
            cout << "ERROR: Unknown spectrum type!" << '\n';
            break;
        }
        }
        break;
    }

    case 2:
    {
        switch (spectrumType)
        {
        case 1:
        {
            const vector<double> dropletSpectrum(spectrumInstance.calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks));
            dropletInstance.setSpectrum(dropletSpectrum);
            dropletInstance.specificStateProperties(startTemperature, spectrumEnergyState, outputPath);
            break;
        }
        case 2:
        {
            unsigned spectrumEnergyStates = 500;
            const vector<double> dropletSpectrum(spectrumInstance.calculateMasterThesisSpectrum(particlesInDroplet, gdd, V, L, spectrumEnergyStates));
            dropletInstance.setSpectrum(dropletSpectrum);
            dropletInstance.specificStateProperties(startTemperature, spectrumEnergyState, outputPath);
            break;
        }
        default:
        {
            cout << "ERROR: Unknown spectrum type!" << '\n';
            break;
        }
        }
        break;
    }

    case 3:
    {
        switch (spectrumType)
        {
        case 1:
        {
            const vector<double> dropletSpectrum(spectrumInstance.calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks));
            dropletInstance.setSpectrum(dropletSpectrum);
            dropletInstance.specificStateTemperatureImpact(spectrumEnergyState, startTemperature, endTemperature, temperatureStep, outputPath);
            break;
        }
        case 2:
        {
            unsigned spectrumEnergyStates = 500;
            const vector<double> dropletSpectrum(spectrumInstance.calculateMasterThesisSpectrum(particlesInDroplet, gdd, V, L, spectrumEnergyStates));
            dropletInstance.setSpectrum(dropletSpectrum);
            dropletInstance.specificStateTemperatureImpact(spectrumEnergyState, startTemperature, endTemperature, temperatureStep, outputPath);
            break;
        }
        default:
        {
            cout << "ERROR: Unknown spectrum type!" << '\n';
            break;
        }
        }
        break;
    }

    default:
    {
        cout << "ERROR: Unknown option!" << '\n';
        break;
    }
    }

    auto tk = high_resolution_clock::now();
    duration<double, std::milli> ms_double = tk - tp;
    cout << "Execution time on the CPU: " << ms_double.count() * 1e-3 << " s";
    outputFile.close();

    return EXIT_SUCCESS;
}
