#include <iostream>
#include <iomanip>
#include <fstream>
#include "duals/dual"
#include "quadrature_multipoly.hpp"
using namespace algoim;
using duals::duald;
using namespace duals::literals;
// Define the function to be differentiated
template <typename T>
duald computeFy(T r_a)
{
    T x_a = 4.0 * r_a * r_a;
    return x_a;
}
/// compute LSF sensitivity
template <int N>
duald computeSens(duald r_a)
{
    // declare active input variable/s
    auto circle = [&](const uvector<duald, 2> &x)
    {
        duald x_a = computeFy(r_a);
        std::cout << "x_a " << x_a << std::endl;
        double yd = 4.0;
        x_a.dpart(yd);
        duald y_a = x_a * x_a + x(1) * x(1) - r_a * r_a;
        return y_a;
    };
    /// general level-set function restricted to 1-D
    auto LSF = [&](const uvector<real, 1> &x)
    {
        uvector<duald, N> xs;
        xs(0) = 2.0;
        xs(1) = x(0);
        duald y = circle(xs);
        return y;
    };
    uvector<real, N - 1> xs;
    xs(0) = 2.0;
    duald y_a = LSF(xs);
    return y_a;
}

int main()
{
    // Evaluate the function at x = 2
    duals::hyperduald x(2 + 0_e, 0 + 1_e);
    // auto y = computeFy(2 + 1_e);
    duald y = computeSens<2>(2 + 1_e);
    // Print the function value and derivative at x = 2
    std::cout << "f(2) = " << y.rpart() << std::endl;
    std::cout << "df/dx(2) = " << y.dpart() << std::endl;
    return 0;
}