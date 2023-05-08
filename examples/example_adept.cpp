#include <iostream>
#include <iomanip>
#include <fstream>
#include "quadrature_multipoly.hpp"
#include "adept.h"
using namespace algoim;
using adept::adouble;

double computeFy(double r_a)
{
    double x_a = 4.0 * r_a * r_a;
    return x_a;
}

/// compute LSF sensitivity
template <int N>
adouble computeSens(adouble r_a)
{
    // declare active input variable/s
    auto circle = [&](const uvector<adouble, 2> &x)
    {
        double xv = computeFy(r_a.value());
        adouble x_a = xv;
        std::cout << "x_a " << x_a << std::endl;
        // stack.new_recording();  // set_gradient works only when I start new recording
        x_a.set_gradient(16.0); // set the gradient of x_a explicitly
        std::cout << "x_a " << x_a << std::endl;
        adouble y_a = x_a * x_a + x(1) * x(1)- r_a * r_a;
        return y_a;
    };
    /// general level-set function restricted to 1-D
    auto LSF = [&](const uvector<real, 1> &x)
    {
        uvector<adouble, N> xs;
        xs(0) = 2.0;
        xs(1) = x(0);
        adouble y = circle(xs);
        return y;
    };
    uvector<real, N - 1> xs;
    xs(0) = 2.0;
    adouble y_a = LSF(xs);
    return y_a;
}

#if ALGOIM_EXAMPLES_DRIVER == 0 || ALGOIM_EXAMPLES_DRIVER == 4
int main(int argc, char *argv[])
{
    // real dwdr = 0.0;
    const int N = 2;
    adouble r_a = 2.0;
    stack.new_recording();
    adouble y_a = computeSens<N>(r_a);
    real dy_dr;
#if 0
    y_a.set_gradient(1.0);
    stack.compute_adjoint();
    r_a.get_gradient(dy_dr);
#endif
#if 1
    r_a.set_gradient(2.0);
    stack.compute_tangent_linear();
    y_a.get_gradient(dy_dr);
#endif
    std::cout << "dy_dr " << dy_dr << std::endl;
    return 0;
}
#endif