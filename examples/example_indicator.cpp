#include <fstream>
#include <iostream>
#include <iomanip> // std::setprecision()
#include "algoim_quad.hpp"
using namespace ::blitz;
using namespace ::Algoim;
bool ellipse = false;
bool circle = false;
bool airfoil = true;
#include <chrono>
using namespace std::chrono;
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
                dX.push_back(0.5 * (Xc.at(i + 1)(j) - Xc.at(i - 1)(j)));
            }
        }
        double nx = dX.at(1);
        double ny = -dX.at(0);
        double ds = 1.0 / sqrt((nx * nx) + (ny * ny));
        nsurf(0) = nx * ds;
        nsurf(1) = ny * ds;
        if (i == nbnd - 1)
        {
            nsurf(0) = 1.0;
            nsurf(1) = 0;
        }
        nor.push_back(nsurf);
    }
    return nor;
}
/// this function constructs the curvature at given boundary points of level-set geometry
template <int N>
std::vector<TinyVector<double, N - 1>> getCurvature(std::vector<TinyVector<double, N>> Xc)
{
    std::vector<TinyVector<double, N - 1>> kappa;
    int nbnd = Xc.size();
    double a0 = 0.2969;
    double a1 = -0.126;
    double a2 = -0.3516;
    double a3 = 0.2843;
    double a4 = -0.1036; // -0.1036 for closed te
    double tc = 0.12;
    double theta = 0.0;
    // thickness
    for (int i = 0; i < nbnd; i++)
    {
        double xc = Xc.at(i)(0);
        double yt = tc * ((a0 * sqrt(xc)) + (a1 * xc) + (a2 * (pow(xc, 2))) + (a3 * (pow(xc, 3))) + (a4 * (pow(xc, 4)))) / (0.2);
        double dytdx = tc * ((0.5 * a0 / sqrt(xc)) + a1 + (2.0 * a2 * xc) + (3.0 * a3 * pow(xc, 2)) + (4.0 * a4 * pow(xc, 3))) / (0.2);
        double d2ytdx = tc * ((-0.25 * a0 / pow(xc, 1.5)) + (2.0 * a2) + (6.0 * a3 * xc) + (12.0 * a4 * pow(xc, 2))) / 0.2;
        double ysu = yt * cos(theta);
        double dydx = dytdx * cos(theta);
        double d2ydx = d2ytdx * cos(theta);
        double roc = (pow((1 + (dydx * dydx)), 1.5)) / abs(d2ydx);
        if (xc == 0.0)
        {
            double rle = 0.5 * pow(a0 * tc / 0.20, 2);
            kappa.push_back(1.0 / rle);
        }
        else if (i == 0 || i == nbnd - 1)
        {
            kappa.push_back(0.0);
        }
        else
        {
            kappa.push_back(1.0 / roc);
        }

        // kappa.push_back(0.0);
    }
    return kappa;
}

template <int N>
struct LevelSet
{
    TinyVector<double, N> norm;
    TinyVector<double, N> xbnd;
    TinyVector<double, N> kappa;
    TinyVector<double, N> tan = {norm(1), -norm(0)};
    template <typename T>
    T operator()(const blitz::TinyVector<T, N> &x) const
    {

        TinyVector<T, N> x_diff;
        x_diff = x - xbnd;
        T dist1 = dot(x_diff, norm);
        double nx = norm(0);
        double ny = norm(1);
        TinyVector<double, N> norvx, norvy;
        norvx = {1.0 - (nx * nx), -nx * ny};
        norvy = {-nx * ny, 1.0 - (ny * ny)};
        TinyVector<T, N> proj;
        proj(0) = dot(norvx, x_diff);
        proj(1) = dot(norvy, x_diff);
        T dist2 = 0.5 * kappa(0) * dot(x_diff, proj);
        T dist = dist1 + dist2;
        return dist;
    }

    template <typename T>
    blitz::TinyVector<T, N> grad(const blitz::TinyVector<T, N> &x) const
    {
        TinyVector<T, N> x_diff;
        x_diff = x - xbnd;
        T perp_tan = dot(tan, x_diff);
        T perp2_x = kappa(0) * tan(0) * perp_tan;
        T perp2_y = kappa(0) * tan(1) * perp_tan;
        return blitz::TinyVector<T, N>(norm(0) + perp2_x, norm(1) + perp2_y);
    }
};
template <int N>
double calcbasis(std::vector<TinyVector<double, N>> xbnd, TinyVector<double, N> xcent,
                 double rho, double delta, double min_dist, int nbnd)
{

    double denom = 0.0;
    /// evaluate the level-set
    for (int idx = 0; idx < nbnd; ++idx)
    {
        TinyVector<double, N> x_diff;
        x_diff = xcent - xbnd.at(idx);
        double delx = sqrt(magsqr(x_diff) + delta);
        double expc = exp(-rho * (delx - min_dist));
        denom += expc;
    }
    return denom;
}

int main(int argc, char *argv[])
{

    if (ellipse)
    {
        const char *area_err = "ellipse_area_err_indicator_lin.dat";
        const char *peri_err = "ellipse_peri_err_indicator_lin.dat";
        ofstream file_area_err, file_peri_err;
        file_area_err.open(area_err, ios::app);
        file_peri_err.open(peri_err, ios::app);
        file_area_err << std::fixed << setprecision(20) << endl;
        file_peri_err << std::fixed << setprecision(20) << endl;
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 2;
        const int count = 6;
        const int N = 2;
        double delta = 1e-10;
        for (int q0 = 0; q0 < count; ++q0)
        {
            /// grid size
            int n = nel[q0];
            /// # boundary points
            int nbnd = n * n;
            cout << "nbnd " << nbnd << endl;
            /// parameters
            double rho = 1 * nbnd;
            cout << "rho " << rho << endl;
            /// major axis
            double a = 4.0;
            /// minor axis
            double b = 1.0;
            double exact_area = M_PI * a * b;
            std::vector<TinyVector<double, N>> Xc;
            std::vector<TinyVector<double, N>> nor;
            std::vector<TinyVector<double, N - 1>> kappa;
            for (int k = 0; k < nbnd; ++k)
            {
                double theta = k * 2.0 * M_PI / nbnd;
                TinyVector<double, N> x, nrm;
                x(0) = (a * cos(theta));
                x(1) = (b * sin(theta));
                // nrm(0) = 2.0 * cos(theta) / a;
                // nrm(1) = 2.0 * sin(theta) / b;
                nrm(0) = 2.0 * (x(0)) / (a * a);
                nrm(1) = 2.0 * (x(1)) / (b * b);
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
                curv(0) = num / den;
                //kappa.push_back(curv);
                kappa.push_back(0.0);
            }
            cout << "qo " << qo1 << endl;
            QuadratureRule<N> qarea, qperi;
            /// Area of a 2D ellipse, computed via the cells of a Cartesian grid
            {
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                /// evaluate levelset and it's gradient
                LevelSet<N> phi;
                TinyVector<double, N> x;
                x(0) = 4.0;
                x(1) = 0.0;
                double dx = 8.2 / n;
                double dy = 2.2 / n;
                double min_x = -4.1;
                double min_y = -1.1;
                double area = 0.0;
                double peri = 0.0;
                double tol = 1e-14;
                double area_weight = 0.0;
                double peri_weight = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        blitz::TinyVector<double, 2> xcent = 0.5 * (xmin + xmax);
                        TinyVector<double, N> dx = xcent - xmin;
                        double L = sqrt(magsqr(dx));

                        // cout << "xmin, xmax " << xmin << "  ,  " << xmax << endl;
                        // cout << "xcent " << xcent << endl;
                        // cout << "L " << L << endl;
                        /// find minimum distance
                        double min_dist = 1e+100;
                        /// find minimum distance
                        for (int idx = 0; idx < nbnd; ++idx)
                        {
                            TinyVector<double, 2> x_diff;
                            x_diff = xcent - Xc.at(idx);
                            double deli_x = sqrt(magsqr(x_diff) + delta);
                            min_dist = min(min_dist, deli_x);
                        }
                        double denom = calcbasis<2>(Xc, xcent, rho, delta, min_dist, nbnd);
                        double fac = exp(2.0 * rho * L);
                        int itr = 0;
                        for (int k = 0; k < nbnd; ++k)
                        {
                            phi.norm = nor.at(k);
                            phi.kappa = kappa.at(k);
                            phi.xbnd = Xc.at(k);
                            TinyVector<double, N> x_diff;
                            x_diff = xcent - Xc.at(k);
                            double delk_x = sqrt(magsqr(x_diff) + delta);
                            double psi = exp(-rho * (delk_x - min_dist)) / denom;
                            if (psi * fac > tol)
                            {
                                itr++;
                                auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                                area += q([](TinyVector<double, 2> x)
                                          { return 1.0; });
                                for (const auto &pt : q.nodes)
                                {
                                    TinyVector<double, N> xp;
                                    xp[0] = pt.x[0];
                                    xp[1] = pt.x[1];
                                    /// find minimum distance
                                    min_dist = 1e+100;
                                    for (int idx = 0; idx < nbnd; ++idx)
                                    {
                                        TinyVector<double, 2> x_diff;
                                        x_diff = xp - Xc.at(idx);
                                        /// this needs to be replaced with interval type
                                        double deli_x = sqrt(magsqr(x_diff) + delta);
                                        min_dist = min(min_dist, deli_x);
                                    }
                                    denom = calcbasis<2>(Xc, xp, rho, delta, min_dist, nbnd);
                                    TinyVector<double, N> x_diff;
                                    x_diff = xp - Xc.at(k);
                                    double delk_x = sqrt(magsqr(x_diff) + delta);
                                    psi = exp(-rho * (delk_x - min_dist)) / denom;
                               
                                    if (psi * pt.w > tol)
                                    {
                                        area_weight += psi * pt.w;
                                        double wt = pt.w * psi;
                                        qarea.evalIntegrand(xp, wt);
                                    }
                                }
                                auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                                for (const auto &pt : qp.nodes)
                                {
                                    TinyVector<double, N> xp;
                                    xp[0] = pt.x[0];
                                    xp[1] = pt.x[1];
                                    /// find minimum distance
                                    min_dist = 1e+100;
                                    for (int idx = 0; idx < nbnd; ++idx)
                                    {
                                        TinyVector<double, 2> x_diff;
                                        x_diff = xp - Xc.at(idx);
                                        double deli_x = sqrt(magsqr(x_diff) + delta);
                                        min_dist = min(min_dist, deli_x);
                                    }
                                    denom = calcbasis<2>(Xc, xp, rho, delta, min_dist, nbnd);
                                    TinyVector<double, N> x_diff;
                                    x_diff = xp - Xc.at(k);
                                    double delk_x = sqrt(magsqr(x_diff) + delta);
                                    psi = exp(-rho * (delk_x - min_dist)) / denom;
                                    
                                    if (psi * pt.w > tol)
                                    {
                                        peri_weight += psi * pt.w;
                                        qperi.evalIntegrand(xp, pt.w);
                                    }
                                }
                            }
                        }
                        // cout << "#local dist functions used: " << itr << endl;
                    }
                cout << "exact_area " << exact_area << endl;
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                std::cout << "  computed area = " << area_weight << "\n";
                double area_err = abs(area_weight - exact_area);
                std::cout << "  area error = " << area_err << "\n";
                std::cout << "  computed perimeter = " << peri_weight << "\n";
                double peri_err = abs(peri_weight - 17.156843550313663);
                std::cout << "  perimeter error = " << peri_err << "\n";
                // member function on the duration object
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                file_area_err << area_err << " ";
                file_peri_err << peri_err << " ";
            }
            std::ofstream farea("ellipse_quad_indicator.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
            std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                      << " quadrature points\n";
            std::ofstream fperi("ellipse_peri_indicator.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
            std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                      << " quadrature points\n";
            // ++qo;
        } /// loop over grid size ends
    }
    if (airfoil)
    {
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 1;
        const int count = 1;
        const int N = 2;
        double delta = 1e-10;
        for (int q0 = 0; q0 < count; ++q0)
        {
            const char *geometry_file = "geometry_data/NACA_0012_200.dat";
            ifstream file;
            file.open(geometry_file);
            std::vector<TinyVector<double, N>> Xcoord;
            std::vector<TinyVector<double, N>> nor;
            std::vector<TinyVector<double, N - 1>> kappa;
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
                Xcoord.push_back(x);
            }
            file.close();
            /// get the number of boundary points
            int nbnd = Xcoord.size();
            cout << "nbnd " << nbnd << endl;
            /// parameters
            double rho = 10 * nbnd;
            cout << "rho " << rho << endl;
            /// construct the normal vector for all boundary points
            nor = constructNormal<N>(Xcoord);
            kappa = getCurvature<N>(Xcoord);
            /// grid size
            int n = nel[1];
            /// translate airfoil
            TinyVector<double, N> xcent;
            xcent(0) = 19.5;
            xcent(1) = 20.0;
            std::vector<TinyVector<double, N>> Xc;
            for (int k = 0; k < nbnd; ++k)
            {
                TinyVector<double, N> xs;
                for (int d = 0; d < N; ++d)
                {
                    xs(d) = Xcoord.at(k)(d) + xcent(d);
                }
                Xc.push_back(xs);
            }
            cout << "qo " << qo1 << endl;
            QuadratureRule<N> qarea, qperi;
            /// Area of a 2D ellipse, computed via the cells of a Cartesian grid
            {
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                /// evaluate levelset and it's gradient
                LevelSet<N> phi;
                TinyVector<double, N> x;
                x(0) = 19.5;
                x(1) = 0.0;
                double dx = 1.2 / n;
                double dy = 0.2 / n;
                double min_x = 19.4;
                double min_y = 19.9;
                double area = 0.0;
                double peri = 0.0;
                double tol = 1e-14;
                double area_weight = 0.0;
                double peri_weight = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        blitz::TinyVector<double, 2> xcent = 0.5 * (xmin + xmax);
                        TinyVector<double, N> dx = xcent - xmin;
                        double L = sqrt(magsqr(dx));

                        // cout << "xmin, xmax " << xmin << "  ,  " << xmax << endl;
                        // cout << "xcent " << xcent << endl;
                        // cout << "L " << L << endl;
                        /// find minimum distance
                        double min_dist = 1e+100;
                        /// find minimum distance
                        for (int idx = 0; idx < nbnd; ++idx)
                        {
                            TinyVector<double, 2> x_diff;
                            x_diff = xcent - Xc.at(idx);
                            double deli_x = sqrt(magsqr(x_diff) + delta);
                            min_dist = min(min_dist, deli_x);
                        }
                        double denom = calcbasis<2>(Xc, xcent, rho, delta, min_dist, nbnd);
                        double fac = exp(2.0 * rho * L);
                        int itr = 0;
                        for (int k = 0; k < nbnd; ++k)
                        {
                            phi.norm = nor.at(k);
                            phi.kappa = kappa.at(k);
                            phi.xbnd = Xc.at(k);
                            TinyVector<double, N> x_diff;
                            x_diff = xcent - Xc.at(k);
                            double delk_x = sqrt(magsqr(x_diff) + delta);
                            double psi = exp(-rho * (delk_x - min_dist)) / denom;
                            if (psi * fac > tol)
                            {
                                itr++;
                                auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                                area += q([](TinyVector<double, 2> x)
                                          { return 1.0; });
                                for (const auto &pt : q.nodes)
                                {
                                    TinyVector<double, N> xp;
                                    xp[0] = pt.x[0];
                                    xp[1] = pt.x[1];
                                    /// find minimum distance
                                    min_dist = 1e+100;
                                    for (int idx = 0; idx < nbnd; ++idx)
                                    {
                                        TinyVector<double, 2> x_diff;
                                        x_diff = xp - Xc.at(idx);
                                        /// this needs to be replaced with interval type
                                        double deli_x = sqrt(magsqr(x_diff) + delta);
                                        min_dist = min(min_dist, deli_x);
                                    }
                                    denom = calcbasis<2>(Xc, xp, rho, delta, min_dist, nbnd);
                                    TinyVector<double, N> x_diff;
                                    x_diff = xp - Xc.at(k);
                                    double delk_x = sqrt(magsqr(x_diff) + delta);
                                    psi = exp(-rho * (delk_x - min_dist)) / denom;
                                    area_weight += psi * pt.w;
                                    if (psi * pt.w > tol)
                                    {
                                        double wt = pt.w * psi;
                                        qarea.evalIntegrand(xp, wt);
                                    }
                                }
                                auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                                for (const auto &pt : qp.nodes)
                                {
                                    TinyVector<double, N> xp;
                                    xp[0] = pt.x[0];
                                    xp[1] = pt.x[1];
                                    /// find minimum distance
                                    min_dist = 1e+100;
                                    for (int idx = 0; idx < nbnd; ++idx)
                                    {
                                        TinyVector<double, 2> x_diff;
                                        x_diff = xp - Xc.at(idx);
                                        double deli_x = sqrt(magsqr(x_diff) + delta);
                                        min_dist = min(min_dist, deli_x);
                                    }
                                    denom = calcbasis<2>(Xc, xp, rho, delta, min_dist, nbnd);
                                    TinyVector<double, N> x_diff;
                                    x_diff = xp - Xc.at(k);
                                    double delk_x = sqrt(magsqr(x_diff) + delta);
                                    psi = exp(-rho * (delk_x - min_dist)) / denom;
                                    peri_weight += psi * pt.w;
                                    if (psi * pt.w > tol)
                                    {
                                        qperi.evalIntegrand(xp, pt.w);
                                    }
                                }
                            }
                        }
                        // cout << "#local dist functions used: " << itr << endl;
                    }
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                std::cout << "  computed area = " << area_weight << "\n";
                std::cout << "  computed perimeter = " << peri_weight << "\n";
                // member function on the duration object
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
            }
            std::ofstream farea("airfoil_area_indicator.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
            std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                      << " quadrature points\n";
            std::ofstream fperi("airfoil_peri_indicator.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
            std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                      << " quadrature points\n";
            // ++qo;
        } /// loop over grid size ends
    }
    return 0;
} // main