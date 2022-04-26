#include <fstream>
#include "algoim_levelset.hpp"
#include <iostream>
using namespace std;

/// this function constructs the normal vectors for given boundary points of level-set geometry
/// \param[in] Xc - boundary coordinates
/// \param[out] nor - boundary normal vectors
template <int N>
std::vector<TinyVector<double, N>> constructNormal(std::vector<TinyVector<double, N>> Xc)
{
    std::vector<TinyVector<double, N>> nor;
    int nbnd = Xc.size();
    TinyVector<double, N> nsurf;
    for (int i = 0; i < nbnd; i++)
    {
        std::vector<double> dX;
        for (int j = 0; j < N; ++j)
        {
            if (i == 0)
            {
                dX.push_back(0.5 * (Xc.at(i + 1)(j) - Xc.at(nbnd - 1)(j)));
            }
            else if (i == nbnd - 1)
            {
                dX.push_back(0.5 * (Xc.at(0)(j) - Xc.at(i - 1)(j)));
            }
            else
            {
                dX.push_back(0.5 * (Xc.at(i + 1)(j) - Xc.at(i - 1))(j));
            }
        }
        double nx = dX.at(1);
        double ny = -dX.at(0);
        double ds = 1.0 / sqrt((nx * nx) + (ny * ny));
        nsurf(0) = nx * ds;
        nsurf(1) = ny * ds;
        nor.push_back(nsurf);
    }
    return nor;
}
/// This functions evaluates the level-set bounds along given x = xslice
/// \param[in] Xc - boundary coordinates
/// \param[in] nor - boundary normal vectors
/// \param[in] kappa - boundary curvature
/// \param[in] xslice - slice in x direction
/// \param[in] ymax - maximum `y` coordinate
/// \param[in] ymin - minimum `y` coordinate
/// \param[in] nel - number of elements in `y` direction
/// \param[in] rho - penalty parameter
/// \param[in] delta- parameter that smooths distance near zero
template <int N>
void testLevelSetBounds(vector<TinyVector<double, N>> Xc, vector<TinyVector<double, N>> nor, vector<TinyVector<double, N-1>> kappa,
                        double xslice, double ymax, double ymin, const int nel, int rho, double delta)
{
    /// create level-set object
    Algoim::LevelSet<N> phi;
    phi.initializeLevelSet(Xc, nor, kappa, rho, 1.0, delta);
    TinyVector<double, N> del;
    TinyVector<Interval<N>, N> phic;
    double ds = (ymax - ymin) / nel;
    del = {0.5 * ds, 0.5 * ds};
    double L = mag(del);
    cout << " ------------------------------------------- " << endl;
    cout << "                       phi_bnd       " << endl;
    cout << " ------------------------------------------- " << endl;
    for (int k = 0; k < nel; ++k)
    {
        TinyVector<double, N> xc;
        TinyVector<double, N> yelem;
        xc(0) = xslice;
        yelem(0) = (ymax - ymin) * (k) / nel + ymin;
        yelem(1) = (ymax - ymin) * (k + 1) / nel + ymin;
        xc(1) = 0.5 * (yelem(0) + yelem(1));
        // xc(0) = 2.0;
        // xc(1) = 0.8;
        Interval<N> x_c = Interval<N>(xc(0));
        Interval<N> y_c = Interval<N>(xc(1));
        x_c.delta() = del;
        y_c.delta() = del;
        phic(0) = x_c;
        phic(1) = y_c;
        phi(phic);
    }
}
/// This functions evaluates the level-set gradient bounds along given x = xslice
/// \param[in] Xc - boundary coordinates
/// \param[in] nor - boundary normal vectors
/// \param[in] xslice - slice in x direction
/// \param[in] ymax - maximum `y` coordinate
/// \param[in] ymin - minimum `y` coordinate
/// \param[in] nel - number of elements in `y` direction
/// \param[in] rho - penalty parameter
/// \param[in] delta- parameter that smooths distance near zero
template <int N>
void testLevelSetGradBounds(vector<TinyVector<double, N>> Xc, vector<TinyVector<double, N>> nor, vector<TinyVector<double, N-1>> kappa,
                            double xslice, double ymax, double ymin, const int nel, int rho, double delta)
{
    /// create level-set object
    Algoim::LevelSet<N> phi;
    phi.initializeLevelSet(Xc, nor, kappa, rho, 1.0, delta);
    TinyVector<double, N> del;
    TinyVector<Interval<N>, N> phic;
    double ds = (ymax - ymin) / nel;
    del = {0.5 * ds, 0.5 * ds};
    double L = mag(del);
    cout << " ------------------------------------------- " << endl;
    cout << "           phi_derivative bnds      " << endl;
    cout << " ------------------------------------------- " << endl;
    for (int k = 0; k < nel; ++k)
    {
        TinyVector<double, N> xc;
        TinyVector<double, N> yelem;
        xc(0) = xslice;
        yelem(0) = (ymax - ymin) * (k) / nel + ymin;
        yelem(1) = (ymax - ymin) * (k + 1) / nel + ymin;
        xc(1) = 0.5 * (yelem(0) + yelem(1));
        Interval<N> x_c = Interval<N>(xc(0));
        Interval<N> y_c = Interval<N>(xc(1));
        x_c.delta() = del;
        y_c.delta() = del;
        phic(0) = x_c;
        phic(1) = y_c;
        phi.grad(phic);
    }
}
/// main code
int main(int argc, char *argv[])
{
    /// dimension
    const int N = 2;
/// if reading from file
#if 0
    //const char *geometry_file = "geometry_data/circle_1.dat";
    const char *geometry_file = "geometry_data/ellipse_1.dat";
    ifstream file;
    file.open(geometry_file);
    /// Vector of normal vectors
    vector<TinyVector<double, N>> nor;
    /// Vector of boundary coordinates
    vector<TinyVector<double, N>> Xc;
    /// read the boundary coordinates from user-provided file
    while (1)
    {
        if (file.eof())
            break;
        TinyVector<double, N> x;
        for (int j = 0; j < N; ++j)
        {
            file >> x(j);
        }
        Xc.push_back(x);
    }
    file.close();
#endif

    /// # boundary poiints
    int nbnd = 32;
    cout << "nbnd " << nbnd << endl;
    /// boundary coordinates
    std::vector<TinyVector<double, N>> Xc;
    /// boundary normal vector coords
    std::vector<TinyVector<double, N>> nor;
    /// boundary curvature
    std::vector<TinyVector<double, N-1>> kappa;
    /// parameters
    double rho = 2.0;
    double delta = 1e-10;
    /// major axis
    double a = 4.0;
    /// minor axis
    double b = 1.0;
    for (int k = 0; k < nbnd; ++k)
    {
        double theta = k * 2.0 * M_PI / nbnd;
        TinyVector<double, N> x, nrm;
        x(0) = a * cos(theta);
        x(1) = b * sin(theta);
        nrm(0) = 2.0 * cos(theta) / a;
        nrm(1) = 2.0 * sin(theta) / b;
        double ds = mag(nrm);
        TinyVector<double, N> ni;
        ni = nrm / ds;
        Xc.push_back(x);
        nor.push_back(ni);
        /// curvature correction
        TinyVector<double, N> dx, d2x;
        dx = {-a * sin(theta), b * cos(theta)};
        d2x = {-a * cos(theta), -b * sin(theta)};
        double num = (dx(0) * d2x(1) - dx(1) * d2x(0));
        double mag_dx = mag(dx);
        double den = mag_dx * mag_dx * mag_dx;
        TinyVector<double, N - 1> curv;
        //curv(0) = num / den;
        curv(0) = 0.0;
        kappa.push_back(curv);
    }
    /// evaluate levelset and it's gradient
    Algoim::LevelSet<N> phi;
    phi.initializeLevelSet(Xc, nor, kappa, rho, 1.0, delta);
    TinyVector<double, N> x;
    x(0) = 2.0;
    x(1) = 0.8;
    std::cout << std::setprecision(10) << std::endl;
    cout << "phi " << phi(x) << endl;
    cout << "grad phi " << phi.grad(x) << endl;
    cout << "hessian phi " << phi.hessian(x) << endl;
    cout << " ------------------------------------------- " << endl;
    /// get the bounds
    const int nel = 20;
    double ymin = -2.0;
    double ymax = 2.0;
    double xslice = 3.0;
    cout << setprecision(12) << endl;
    testLevelSetBounds<N>(Xc, nor, kappa, xslice, ymax, ymin, nel, rho, delta);
    //testLevelSetGradBounds<N>(Xc, nor, kappa, xslice, ymax, ymin, nel, rho, delta);
    return 0;
} // main ends
