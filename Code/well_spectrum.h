#ifndef WELL_SPECTRUM_H
#define WELL_SPECTRUM_H
#include <vector>

// This class is used to calculate energy spectrum in finite well potential problem
class wellspec
{
private:
    // Accuracy of the bisection method
    double eps;
    int steps;

    // Pointer to function
    typedef double (wellspec::*function)(double x, double L, double a, double cons);

    // Function used to calculate roots of function to even bond energy levels
    double even_eq_bond(double x, double L, double a, double cons);
    // Function used to calculate roots of function to odd bond energy levels
    double odd_eq_bond(double x, double L, double a, double cons);
    // Function used to calculate roots of function to even scatter energy levels
    double even_eq_scatt(double x, double L, double a, double cons);
    // Function used to calculate roots of function to odd scatter energy levels
    double odd_eq_scatt(double x, double L, double a, double cons);

    // Function used to calculate one root of function
    double find_root(double xp, double xk, function fptr, double L, double a, double cons);
    // Function used to calculate several roots of function
    std::vector<double> find_all_roots(double xp, double end, function fptr, double L, double a, double cons);

public:
    // Constructor, it sets proper accuracy of bisection method and number of steps
    wellspec(double this_eps, int this_steps);
    // Main function used to calculate whole spectrum of finite potential problem
    std::vector<double> spectrum(double L, double a, double V, double cons, double xp, double xkb, double xks);
    // Function to set up number of steps
    void set_steps(int this_steps);
};

#endif // WELL_SPECTRUM_H