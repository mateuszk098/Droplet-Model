/** Author: Mateusz Kowalczyk
 * This program calculates interesting properties of a new quantum state,
 * quantum droplet. The main method of calculating is based on performing
 * the calculations of canonical ensemble based on grand canonical ensemble.
 * For this purpose, I use the path integration formula based on Cauchy's 
 * integration formula. This program using multithreads.
 **/

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include "droplet.h"
#include "well_spectrum.h"
#include "real_spectrum.h"

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(0);
    cin.tie(0);

    // ----------------------- SECTION TO SET UP PROPER PARAMETERS --------------------------
    // Droplet parameters
    double gdd = 25.0;                             // gdd coefficient
    int N_tot = 50;                                // Number of bosons
    double N_droplet = static_cast<double>(N_tot); // Number of bosons whose create droplet

    // Model parameters
    double a = 2.0 * M_PI * M_PI * N_tot / (3.0 * gdd); // Initial droplet width
    double V = 3.0 / (8.0 * M_PI * M_PI) * gdd * gdd;   // Potential
    double L = 10.0 * a;                                // Walls
    double c = 2.0 * V * a * a;                         // Constant coefficient in model

    // Parameters used in finite potential well
    int steps = 10000;          // Number of sub-bisection steps
    double eps = 1e-6;          // Accuracy of sub-bisection method
    double xp = eps;            // Begin of range
    double xkb = sqrt(c) - eps; // End of range for bound levels
    double xks = 50.0 - eps;    // End of range for scatter levels

    // Initialise objects
    vector<double> full_spectrum;            // Vector to store energy values
    droplet *d = new droplet(N_tot, a, V);   // Create new initial droplet
    wellspec *ws = new wellspec(eps, steps); // Create new spectrum from finite well
    realspec *rs = new realspec;
    // ----------------------- END OF SECTION TO SET UP PROPER PARAMETERS -------------------

    // ----------------------- SIMULATION LOOP ----------------------------------------------

    string filename = "../Data/";
    filename.append(argv[1]);
    ofstream output_file(filename.c_str(), ios::out);
    int thermalization = 5;                 // Self-consistent calculation steps
    auto tp = high_resolution_clock::now(); // Execution time - start

    // Main temperature loop
    for (double T = 0.1; T < 30.0; T += 0.1) // Iterate by temperature
    {
        double beta = 1.0 / T; // Beta coefficient

        // Thermalization loop
        // N_drop changes by reference 5 times
        // If thermalization reach 5th step then output about droplet is saved
        for (int index = 0; index < thermalization; index++)
        {
            // -------------- POTENTIAL WELL SPECTRUM ---------------------------------------
            // Parameters a, c and xkb change with every step
            a = 2.0 * M_PI * M_PI * N_droplet / (3.0 * gdd);
            c = 2.0 * V * a * a;
            xkb = sqrt(c) - eps;
            full_spectrum = ws->spectrum(L, a, V, c, xp, xkb, xks);
            // -------------- END OF POTENTIAL WELL SPECTRUM --------------------------------

            // -------------- REAL SPECTRUM FROM MASTER THESIS ------------------------------
            // double scatter_levels = 1000; // Number of scatter levels from real spectrum
            // full_spectrum = rs->spectrum(N_droplet, gdd, V, L, scatter_levels);
            // -------------- END OF REAL SPECTRUM FROM MASTER THESIS -----------------------

            d->set_spectrum(full_spectrum);
            string output = (d->droplet_width(beta, N_droplet)).str();

            if (index == thermalization - 1)
                output_file << output;
        }
    }

    output_file.close();
    auto tk = high_resolution_clock::now();           // Execution time - stop
    duration<double, std::milli> ms_double = tk - tp; // Getting  milliseconds as a double
    cout << "Execution time on the CPU: " << ms_double.count() << " ms";
    // ----------------------- END OF SIMULATION LOOP ---------------------------------------

    // ----------------------- TEST OTHER FUNCTIONS -----------------------------------------
    // full_spectrum = ws->spectrum(L, a, V, c, xp, xkb, xks);   // Quantum well
    // full_spectrum = rs->spectrum(N_droplet, gdd, V, L, 1000); // Real
    // d->set_spectrum(full_spectrum);
    // d->specific_state_properties(100, 0, "out.txt");
    // d->state_with_temperature_change(0, 1.0, 100.0, 5.0, "out.txt");
    // ----------------------- END OF TEST OTHER FUNCTIONS ----------------------------------

    delete rs;
    delete ws;
    delete d;

    return 0;
}
