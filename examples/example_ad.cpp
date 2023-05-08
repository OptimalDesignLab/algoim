#include <iostream>
#include <iomanip>
#include <fstream>

#include "adept.h"
using adept::adouble;
adept::Stack stack;

double computeFy(double r_a)
{
    double x_a = 4.0 * r_a * r_a;
    return x_a;
}

adouble computeSens(adouble r_a)
{
    double x_a_passive = computeFy(r_a.value());
    adouble x_a;
    // Assign to x_a without adding derivative to the stack
    x_a.set_value(x_a_passive);
    // Set the gradient of x_a explicitly
    x_a.add_derivative_dependence(r_a, 4.0);
    adouble y_a = x_a * x_a + 4.0 - r_a * r_a;
    return y_a;
}

int main(int argc, char *argv[])
{
    adouble r_a = 2.0;
    stack.new_recording();
    adouble y_a = computeSens(r_a);
    double dy_dr;
    r_a.set_gradient(1.0);
    stack.compute_tangent_linear();
    y_a.get_gradient(dy_dr);
    std::cout << "dy_dr " << dy_dr << std::endl;
    stack.print_gradients();
    return 0;
}

