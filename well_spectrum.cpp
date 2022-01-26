#include <algorithm>
#include <cmath>
#include "well_spectrum.h"
using namespace std;

/** Constructor, it sets proper accuracy of bisection method and number of steps
 * @param double eps (accuracy of sub-bisection method)
 * @param int steps (number of sub-bisections in whole range (xp; xk))
 * **/
wellspec::wellspec(double this_eps, int this_steps)
{
    eps = this_eps;
    steps = this_steps;
}

/** This function is used to set number of sub-bestions
 * @param int steps (number of sub-bisections in whole range (xp; xk))
 * **/
void wellspec::set_steps(int this_steps)
{
    steps = this_steps;
}

/** This method returns function to roots of even bound solutions
 * @param double x (function parameter)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of even bound solutions)
 * **/
double wellspec::even_eq_bond(double x, double L, double a, double c)
{
    return x / tanh((1 - (L / a)) * x) + sqrt(c - x * x) * tan(sqrt(c - x * x));
}

/** This method returns function to roots of odd bound solutions
 * @param double x (function parameter)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of odd bound solutions)
 * **/
double wellspec::odd_eq_bond(double x, double L, double a, double c)
{
    return x / tanh((1 - (L / a)) * x) - sqrt(c - x * x) / tan(sqrt(c - x * x));
}

/** This method returns function to roots of even scatter solutions
 * @param double x (function parameter)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of even scatter solutions)
 * **/
double wellspec::even_eq_scatt(double x, double L, double a, double c)
{
    return x / tan((1 - (L / a)) * x) + sqrt(c + x * x) * tan(sqrt(c + x * x));
}

/** This method returns function to roots of odd scatter solutions
 * @param double x (function parameter)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double (function to roots of odd scatter solutions)
 * **/
double wellspec::odd_eq_scatt(double x, double L, double a, double c)
{
    return x / tan((1 - (L / a)) * x) - sqrt(c + x * x) / tan(sqrt(c + x * x));
}

/** This method using bisection method to find root of function in specific subrange (xp; xk)
 * @param double xp (begin of the range)
 * @param double xk (end of the range) 
 * @param function fptr (pointer to appropriate function)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return double root (root of function in the range (xp; xk), 
 *  if method not found root in the range, then -1 is returned)
 * **/
double wellspec::find_root(double xp, double xk, function fptr, double L, double a, double c)
{
    double root = 0; // Searched root of function defined as fptr

    // Dereferences a member function pointer - the parens around this->*fptr are important for precedence.
    // (this->*fptr)(xp, L, a, c);

    // We are looking for the root in range (xp, xk)
    // Bisection method assumes root not exist if f(tp) * f(tk) > 0
    if ((this->*fptr)(xp, L, a, c) * (this->*fptr)(xk, L, a, c) > 0)
        return -1;
    // Else possibly we have root there
    else
    {
        while ((xk - xp) * 0.5 > eps)
        {
            root = (xp + xk) * 0.5;

            if ((this->*fptr)(root, L, a, c) == 0)
                break;
            else if ((this->*fptr)(xp, L, a, c) * (this->*fptr)(root, L, a, c) < 0)
                xk = root;
            else
                xp = root;
        }
    }

    return root;
}

/** This method using multi bisection method to find every roots of function in whole range (xp; xk)
 * @param double xp (begin of the range)
 * @param double xk (end of the range) 
 * @param function fptr (pointer to appropriate function)
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @return std::vector<double> roots (roots of function in the whole range (xp; xk)
 * **/
vector<double> wellspec::find_all_roots(double xp, double xk, function fptr, double L, double a, double c)
{
    double interval = abs(xk - xp) / steps;
    vector<double> roots;
    double root = 0;

    for (int i = 0; i < steps; i++)
    {
        root = find_root(xp, xp + interval, fptr, L, a, c);
        // Omit -1 roots might be weird but we know we are looking
        // for only positive roots so this is a correct approach
        // Apart from that, function find_root must returns something
        // So I adopt if root does not exist, -1 will be returned
        // But the interval must be increased
        if (root == -1)
        {
            xp += interval;
            continue;
        }
        roots.push_back(root);
        xp += interval;
    }

    return roots;
}

/** This method is used to find solutions of finite well potential problem
 * @param double L (half-width of the infinite wall)
 * @param double a (half-width of the finite well)
 * @param double V (depth of the finite well)
 * @param double c (constant parameter characterising the well c = (2m/hbar^2)V*a^2)
 * @param double xp (begin of the range)
 * @param double xk (end of the range) 
 * @return std::vector<double> solutions of finite well potential problem
 * **/
vector<double> wellspec::spectrum(double L, double a, double V, double c, double xp, double xkb, double xks)
{
    vector<double> even_bond_roots;
    vector<double> odd_bond_roots;
    vector<double> even_scatt_roots;
    vector<double> odd_scatt_roots;
    vector<double> bond_roots;
    vector<double> scatt_roots;
    vector<double> full_spectrum;

    // Calculate roots of even bound function
    even_bond_roots = find_all_roots(xp, xkb, even_eq_bond, L, a, c);
    // Calculate roots of even odd function
    odd_bond_roots = find_all_roots(xp, xkb, odd_eq_bond, L, a, c);

    // Now I have to change range of solutions because scatter levels do not have
    // end of range condition
    // Calculate roots of even scatter function
    even_scatt_roots = find_all_roots(xp, xks, even_eq_scatt, L, a, c);
    // Calculate roots of odd scatter function
    odd_scatt_roots = find_all_roots(xp, xks, odd_eq_scatt, L, a, c);

    // I have to remove roots that have been incorrectly identified as roots, because they are asymptotes
    // I just check if function in this root is greater than 1, if so, it is rather the asymptote
    // I will use appropriate lambda functions

    // Asymptotes in even bound equation
    auto lambda_ebr = [this, L, a, c](double x)
    { return abs(even_eq_bond(x, L, a, c)) > 1.0; };
    // Asymptotes in odd bound equation
    auto lambda_obr = [this, L, a, c](double x)
    { return abs(odd_eq_bond(x, L, a, c)) > 1.0; };
    // Asymptotes in even scatter equation
    auto lambda_esr = [this, L, a, c](double x)
    { return abs(even_eq_scatt(x, L, a, c)) > 1.0; };
    // Asymptotes in odd scatter equation
    auto lambda_osr = [this, L, a, c](double x)
    { return abs(odd_eq_scatt(x, L, a, c)) > 1.0; };

    // Now I remove these asymptotes from appropriate vectors
    even_bond_roots.erase(remove_if(even_bond_roots.begin(), even_bond_roots.end(), lambda_ebr),
                          even_bond_roots.end());
    odd_bond_roots.erase(remove_if(odd_bond_roots.begin(), odd_bond_roots.end(), lambda_obr),
                         odd_bond_roots.end());
    even_scatt_roots.erase(remove_if(even_scatt_roots.begin(), even_scatt_roots.end(), lambda_esr),
                           even_scatt_roots.end());
    odd_scatt_roots.erase(remove_if(odd_scatt_roots.begin(), odd_scatt_roots.end(), lambda_osr),
                          odd_scatt_roots.end());

    // I combine bound solutions
    bond_roots.insert(bond_roots.begin(), even_bond_roots.begin(), even_bond_roots.end());
    bond_roots.insert(bond_roots.end(), odd_bond_roots.begin(), odd_bond_roots.end());
    // I combine scatter solutions
    scatt_roots.insert(scatt_roots.begin(), even_scatt_roots.begin(), even_scatt_roots.end());
    scatt_roots.insert(scatt_roots.end(), odd_scatt_roots.begin(), odd_scatt_roots.end());

    // I evaluate bound energy levels
    auto lambda_bound_coeff = [a](double &x)
    { return -(x * x) / (2 * a * a); };
    transform(bond_roots.begin(), bond_roots.end(), bond_roots.begin(), lambda_bound_coeff);

    // I evaluate scatter energy levels
    auto lambda_scatter_coeff = [a](double &x)
    { return (x * x) / (2 * a * a); };
    transform(scatt_roots.begin(), scatt_roots.end(), scatt_roots.begin(), lambda_scatter_coeff);

    // I combine bound and scatter spectrum and sort
    full_spectrum.insert(full_spectrum.begin(), bond_roots.begin(), bond_roots.end());
    full_spectrum.insert(full_spectrum.end(), scatt_roots.begin(), scatt_roots.end());
    sort(full_spectrum.begin(), full_spectrum.end());

    // Eventually I translate a whole spectrum by potential depth
    auto lambda_translate = [V](double &x)
    { return x + V; };
    transform(full_spectrum.begin(), full_spectrum.end(), full_spectrum.begin(), lambda_translate);

    // And I add first droplet energy state equal zero
    full_spectrum.insert(full_spectrum.begin(), 0);

    return full_spectrum;
}
