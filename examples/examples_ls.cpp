#include <fstream>
#include "/users/kaurs3/Sharan/Research/Spring_2020/algoim_fork/src/algoim_levelset.hpp"
#include <iostream>
// #include "algoim_quad.hpp"
int ls = 2;
/// this function constructs the normal vectors for given boundary points of level-set geometry
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
                dX.push_back(Xc.at(i + 1)(j) - Xc.at(i)(j));
            }
            else if (i == nbnd - 1)
            {
                dX.push_back(Xc.at(i)(j) - Xc.at(i - 1)(j));
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
template <int N>
struct Ellipsoid
{
    template <typename T>
    T operator()(const blitz::TinyVector<T, N> &x) const
    {
        return x(0) * x(0) + 16.0 * x(1) * x(1) - 16.0;
    }

    template <typename T>
    blitz::TinyVector<T, N> grad(const blitz::TinyVector<T, N> &x) const
    {
        return blitz::TinyVector<T, N>(2.0 * x(0), 32.0 * x(1));
    }
};

template <int N>
struct Circle
{
    template <typename T>
    T operator()(const blitz::TinyVector<T, N> &x) const
    {
        return ((x(0)-2.0) * (x(0) - 2.0)) + ((x(1)-2.0)* (x(1)-2.0)) - 1.0;
    }

    template <typename T>
    blitz::TinyVector<T, N> grad(const blitz::TinyVector<T, N> &x) const
    {
        return blitz::TinyVector<T, N>(2.0 * (x(0)-2.0), 2.0*(x(1)-2.0));
    }
};

int main(int argc, char *argv[])
{
    /// dimension
    const int N = 2;
    //const char *geometry_file = "geometry_data/circle_1.dat";
    const char *geometry_file = "geometry_data/ellipse.dat";
    ifstream file;
    file.open(geometry_file);
    /// Vector of normal vectors
    std::vector<TinyVector<double, N>> nor;
    /// Vector of boundary coordinates
    std::vector<TinyVector<double, N>> Xc;
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
    /// get the number of boundary points
    int nbnd = Xc.size();
    cout << "nbnd " << nbnd << endl;
    /// construct the normal vector for all boundary points
    nor = constructNormal<N>(Xc);
    /// parameters
    double rho = 10 * nbnd;
    double delta = 1e-10;
    /// evaluate levelset and it's gradient
    Algoim::LevelSet<N> phi(Xc, nor, rho, delta);
    TinyVector<double, N> x;
    x(0) = 2;
    x(1) = 0.844;
    std::cout << std::setprecision(10) << std::endl;
    cout << "phi " << phi(x) << endl;
    cout << "grad phi " << phi.grad(x) << endl;
    cout << "hessian phi " << phi.hessian(x) << endl;
    cout << " ------------------------------------------- " << endl;


    /// Area of a 2D ellipse, computed via the cells of a Cartesian grid
    // {
    //     int n = 64;
    //     std::cout << "Area of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
    //     double dx = 8.2 / n;
    //     double dy = 2.2 / n;
    //     double min_x = -4.1;
    //     double min_y = -1.1;
    //     double area = 0.0;
    //     for (int i = 0; i < n; ++i)
    //         for (int j = 0; j < n; ++j)
    //         {
    //             blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
    //             blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
    //             //cout << "xmin " << xmin  << "     " << "xmax " << xmax << endl;
    //             auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, 4);
    //             area += q([](TinyVector<double, 2> x) { return 1.0; });
    //         }
    //     std::cout << "  computed area = " << area << "\n";
    //     std::cout << "    exact area = " << 4.0 * M_PI << "\n";
    // }
    #if 0
    //"Area of a 2D circle using automatic subdivision
    {
        std::cout << "Area of a 2D circle using automatic subdivision:\n";
        blitz::TinyVector<double, N> xupper;
        blitz::TinyVector<double, N> xlower;
        xlower = {0.9, 0.9};
        xupper = {3.1, 3.1};
        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, 2);
        double area = q([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed area = " << area << "\n";
        std::cout << "    (exact area = " <<  M_PI << "\n";
        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, 2);
        double peri = qp([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed perimeter = " << peri << "\n";
        std::cout << "    (exact perimeter = " <<  2.0*M_PI << "\n";
    }
    //"Area of a 2D ellipse using automatic subdivision
    {
        Circle<2> phi2;
        std::cout << "Area of a 2D circle using (exact level-set) automatic subdivision:\n";
        blitz::TinyVector<double, N> xupper;
        blitz::TinyVector<double, N> xlower;
        xlower = {0.9, 0.9};
        xupper = {3.1, 3.1};
        auto q = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, 2);
        double area = q([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed area = " << area << "\n";
        std::cout << "    (exact area = " <<  M_PI << "\n";
        auto qp = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, 2);
        double peri = qp([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed perimeter = " << peri << "\n";
        std::cout << "    (exact perimeter = " <<  2.0*M_PI << "\n";
    }
    #endif

    //"Area of a 2D ellipse using automatic subdivision
    {
        std::cout << "Area of a 2D ellipse using automatic subdivision:\n";
        blitz::TinyVector<double, N> xupper;
        blitz::TinyVector<double, N> xlower;
        xlower = {-4.1, -1.1};
        xupper = {4.1, 1.1};
        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, 4);
        double area = q([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed area = " << area << "\n";
        std::cout << "    exact area = " << 4.0 * M_PI << "\n";
        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, 4);
        double peri = qp([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed perimeter = " << peri << "\n";
        std::cout << "    exact perimeter = " <<  17.15684 << "\n";
    }

    //"Area of a 2D ellipse using automatic subdivision
    {
        Ellipsoid<2> phi2;
        std::cout << "Area of a 2D ellipse using (exact level-set) automatic subdivision:\n";
        blitz::TinyVector<double, N> xupper;
        blitz::TinyVector<double, N> xlower;
        xlower = {-4.1, -1.1};
        xupper = {4.1, 1.1};
        auto q = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, 4);
        double area = q([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed area = " << area << "\n";
        std::cout << "    exact area = " << 4.0 * M_PI << "\n";
        auto qp = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, 4);
        double peri = qp([](const TinyVector<double, 2> x) { return 1.0; });
        std::cout << "  computed perimeter = " << peri << "\n";
        std::cout << "    exact perimeter = " <<  17.15684 << "\n";
    }

    return 0;
} // main