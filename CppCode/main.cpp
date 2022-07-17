/** Author: Mateusz Kowalczyk
 * This program calculates interesting properties of a new quantum state,
 * quantum droplet. The main method of calculating is based on performing
 * the calculations of canonical ensemble based on grand canonical ensemble.
 * For this purpose, I use the path integration formula based on Cauchy's
 * integration formula. This program using multithreads.
 **/

#define _USE_MATH_DEFINES

#include "droplet.h"
#include "Spectrum.h"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

void load_parameters(string &input_name, string &output_name, int &spec_type, double &gdd, int &N_tot, double &Tp, double &Tk, double &dT, int &level)
{
    string line, tmp;
    ifstream input;
    input.open(input_name, ios::in);
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> tmp >> output_name >> tmp >> spec_type >> tmp >> gdd >> tmp >> N_tot;
            input >> tmp >> Tp >> tmp >> Tk >> tmp >> dT >> tmp >> level;
        }
        input.close();
    }
    else
        cerr << "Cannot open file!" << '\n';
}

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(0);
    cin.tie(0);

    if (argc < 2)
    {
        cerr << "Too few arguments!" << '\n';
        return EXIT_FAILURE;
    }

    // ----------------------- SECTION TO SET UP PROPER PARAMETERS --------------------------
    // Droplet parameters
    double gdd = 0; // gdd coefficient
    int N_tot = 0;  // Number of bosons

    // Temperature parameters;
    double Tp = 0; // Initial temperature
    double Tk = 0; // Final temperature
    double dT = 0; // Interval

    // Energy level
    int level = 0;
    // Spectrum type
    int spec_type = 0;

    // Filenames
    string input_name = argv[1]; // File with parameters
    string output_name = "";     // File with output

    // Load parameters from file;
    load_parameters(input_name, output_name, spec_type, gdd, N_tot, Tp, Tk, dT, level);

    double N_droplet = static_cast<double>(N_tot); // Number of bosons whose create droplet

    // Model parameters
    double a = 2.0 * M_PI * M_PI * N_tot / (3.0 * gdd); // Initial droplet width
    double V = 3.0 / (8.0 * M_PI * M_PI) * gdd * gdd;   // Potential
    double L = 10. * a;                                 // Walls
    double c = 2.0 * V * a * a;                         // Constant coefficient in model

    // Parameters used in finite potential well
    int steps = 20000;                 // Number of sub-bisection steps
    double eps = 1e-6;                 // Accuracy of sub-bisection method
    double xp = eps;                   // Begin of range
    double xkb = sqrt(c * 0.25) - eps; // End of range for bound levels
    double xks = 100.0 - eps;          // End of range for scatter levels

    // Initialise objects
    vector<double> full_spectrum;            // Vector to store energy values
    Droplet *d = new Droplet(N_tot, a, V);   // Create new initial droplet
    Spectrum *ws = new Spectrum(eps, steps); // Create new spectrum from finite well

    // ----------------------- END OF SECTION TO SET UP PROPER PARAMETERS -------------------

    // ----------------------- SIMULATION LOOP ----------------------------------------------
    string path = "../Data/";
    path.append(output_name);
    ofstream output_file(path.c_str(), ios::out);
    output_file.precision(5);
    int thermalization = 5;                 // Self-consistent calculation steps
    auto tp = high_resolution_clock::now(); // Execution time - start

    // Main temperature loop
    // for (double T = Tp; T < Tk; T += dT) // Iterate by temperature
    // {
    //     // Thermalization loop
    //     // N_drop changes by reference 5 times
    //     // If thermalization reach 5th step then output about droplet is saved
    //     for (int index = 0; index < thermalization; index++)
    //     {
    //         // -------------- POTENTIAL WELL SPECTRUM ---------------------------------------
    //         // Parameters a, c and xkb change with every step
    //         if (spec_type == 1)
    //         {
    //             a = 2.0 * M_PI * M_PI * N_droplet / (3.0 * gdd);
    //             c = 2.0 * V * a * a;
    //             xkb = sqrt(c * 0.25) - eps;
    //             full_spectrum = ws->calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks);
    //         }
    //         // -------------- END OF POTENTIAL WELL SPECTRUM --------------------------------

    //         // -------------- REAL SPECTRUM FROM MASTER THESIS ------------------------------
    //         else if (spec_type == 2)
    //         {
    //             int scatter_levels = 1000; // Number of scatter levels from real spectrum
    //             full_spectrum = ws->calculateMasterThesisSpectrum(N_droplet, gdd, V, L, scatter_levels);
    //         }
    //         // -------------- END OF REAL SPECTRUM FROM MASTER THESIS -----------------------

    //         d->setSpectrum(full_spectrum);
    //         string output = (d->dropletWidth(T, N_droplet)).str();

    //         if (index == thermalization - 1)
    //             output_file << output;
    //     }
    // }

    output_file.close();
    auto tk = high_resolution_clock::now();           // Execution time - stop
    duration<double, std::milli> ms_double = tk - tp; // Getting  milliseconds as a double
    cout << "Execution time on the CPU: " << ms_double.count() * 1e-3 << " s";
    // ----------------------- END OF SIMULATION LOOP ---------------------------------------

    // ----------------------- TEST OTHER FUNCTIONS -----------------------------------------
    // xkb = sqrt(c * 0.25) - eps;
    full_spectrum = ws->calculateWellSpectrum(L * 0.5, a * 0.5, V, c * 0.25, xp, xkb, xks); // Quantum well
    // full_spectrum = rs->spectrum(N_droplet, gdd, V, L, 1000); // Real
    // for (int i = 0; i < 1000; i++) // 1D condensate spectrum
    //     full_spectrum.push_back(i);

    d->setSpectrum(full_spectrum);
    d->specificStateProperties(Tp, level, output_name);
    d->specificStateTemperatureImpact(level, Tp, Tk, dT, output_name);
    // ----------------------- END OF TEST OTHER FUNCTIONS ----------------------------------

    delete ws;
    delete d;

    return EXIT_SUCCESS;
}
