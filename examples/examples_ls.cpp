#include <fstream>
#include "/users/kaurs3/Sharan/Research/Spring_2020/algoim_fork/src/algoim_levelset.hpp"
#include <iostream>

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

int main(int argc, char *argv[])
{
    /// dimension
    const int N = 2;
    const char *geometry_file = "geometry_data/circle_1.dat";
    ifstream file;
    file.open(geometry_file);
    /// Vector of normal vectors
    std::vector<TinyVector<double, N>> nor;
    /// Vector of boundary coordinates
    std::vector<TinyVector<double, N>> Xc;
    /// read the boundary coordinates from user-provided file
    while (1)
    {
        TinyVector<double, N> x;
        for (int j = 0; j < N; ++j)
        {
            file >> x(j);
        }
        Xc.push_back(x);
        if (file.eof())
            break;
    }
    file.close();
    /// get the number of boundary points
    int nbnd = Xc.size();
    /// construct the normal vector for all boundary points
    nor = constructNormal<N>(Xc);
    /// parameters
    double rho = 10 * nbnd;
    double delta = 1e-10;
    /// evaluate levelset and it's gradient
    Algoim::LevelSet<N> phi(Xc, nor, rho, delta);
    TinyVector<double, N> x;
    x(0) = 4;
    x(1) = 0;
    cout << "phi " << phi(x) << endl;
    cout << "grad phi " << phi.grad(x) << endl;
    return 0;
} // main