#include <fstream>
#include "algoim_levelset.hpp"
#include <iostream>
#include <iomanip> // std::setprecision()
bool ellipse = false;
bool circle = false;
bool airfoil = false;
bool airfoil_up = true;
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
        // if (i == 0 || i == nbnd - 1)
        // {
        //     nsurf(0) = 1.0;
        //     nsurf(1) = 0;
        // }
        if (i == 0)
        {
            nsurf(0) = 1.0;
            nsurf(1) = 0;
        }
        nor.push_back(nsurf);
    }
    return nor;
}
/// this function constructs the curvature at given boundary points of level-set geometry
template <int N, bool cte = false>
std::vector<TinyVector<double, N - 1>> getCurvature(std::vector<TinyVector<double, N>> Xc)
{
    std::vector<TinyVector<double, N - 1>> kappa;
    int nbnd = Xc.size();
    double a0 = 0.2969;
    double a1 = -0.126;
    double a2 = -0.3516;
    double a3 = 0.2843;
    double a4 = -0.1015;
    if (cte)
    {
        a4 = -0.1036;
    }
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
        // else if (i == 0 || i == 1 || i == nbnd - 1 || i == nbnd - 2)
        // {
        //     kappa.push_back(0.0);
        // }
        else if (i == 0 || i == 1)
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

int main(int argc, char *argv[])
{
    /// dimension
    const int N = 2;
    double delta = 1e-10;
    if (airfoil)
    {
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 10;
        const int count = 1;
        double exact_area = 0.0817073;
        const int npts = 200;
        const int nbnd = 2 * npts - 1;
        for (int q0 = 0; q0 < count; ++q0)
        {
            std::vector<TinyVector<double, N>> Xc;
            std::vector<TinyVector<double, N>> nor;
            std::vector<TinyVector<double, N - 1>> kappa;
            /// calculate boundary coordinates
            double tc = 0.12;
            TinyVector<double, npts - 1> beta;
            double beta_max = M_PI / 1.022;
            double dbeta = beta_max / (npts - 2);
            beta(0) = beta_max;
            for (int i = 1; i < npts - 1; ++i)
            {
                beta(i) = beta(i - 1) - dbeta;
            }
            TinyVector<double, nbnd> xb;
            TinyVector<double, nbnd> yb;
            constexpr bool cte = true;
            double a0 = 0.2969;
            double a1 = -0.126;
            double a2 = -0.3516;
            double a3 = 0.2843;
            double a4 = -0.1015;
            if (cte)
            {
                a4 = -0.1036;
            }
            /// upper boundary
            for (int i = 0; i < npts; ++i)
            {
                if (i == 0)
                {
                    xb(i) = 1.0;
                }
                else
                {
                    xb(i) = (1.0 - cos(beta(i - 1))) / 2.0;
                }
                double term1 = (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
                double term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
                yb(i) = 5.0 * tc * (term1 + term2);
            }
            /// lower boundary
            for (int i = 0; i < npts - 1; ++i)
            {
                xb(i + npts) = xb(npts - 2 - i);
                yb(i + npts) = -yb(npts - 2 - i);
            }
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x;
                x(0) = xb(i);
                x(1) = yb(i);
                Xc.push_back(x);
            }
            /// get the number of boundary points
            // nbnd = Xc.size();
            cout << "nbnd " << nbnd << endl;
            cout << "Xc.size() " << Xc.size() << endl;
            /// parameters
            int ratio = 10;
            double rho = ratio * nbnd;
            cout << "rho " << rho << endl;
            /// construct the normal vector for all boundary points
            nor = constructNormal<N>(Xc);
            kappa = getCurvature<N, cte>(Xc);
            // cout << "Airfoil coordinates: " << endl;
            // for (int k = 0; k < nbnd; ++k)
            // {
            //     cout << Xc.at(k) << endl;
            // }
            /// grid size
            int n = nel[0];
            /// translate airfoil
            TinyVector<double, N> xcent;
            xcent(0) = 19.5;
            xcent(1) = 20.0;
            std::vector<TinyVector<double, N>> Xcoord;
            for (int k = 0; k < nbnd; ++k)
            {
                TinyVector<double, N> xs;
                for (int d = 0; d < N; ++d)
                {
                    xs(d) = Xc.at(k)(d) + xcent(d);
                }
                Xcoord.push_back(xs);
            }
            /// evaluate levelset and it's gradient
            Algoim::LevelSet<N> phi;
            phi.initializeLevelSet(Xcoord, nor, kappa, rho, 1.0, delta);
            TinyVector<double, N> x;
            x(0) = 19.5;
            x(1) = 20.0;
            cout << "phi " << phi(x) << endl;
            cout << "grad phi " << phi.grad(x) << endl;
            cout << "hessian phi " << phi.hessian(x) << endl;
            cout << " ------------------------------------------- " << endl;
            cout << "qo " << qo1 << endl;
            QuadratureRule<N> qarea, qperi, qarea_auto, qperi_auto;
#if 1
            /// Area and Perimeter of a 2D airfoil, computed via the cells of a Cartesian grid
            {
                const char *surface_quad_file = "quad_points_peri_auto.dat";
                ofstream file;
                file.open(surface_quad_file);
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D airfoil, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double dx = 1.1 / n;
                double dy = 0.2 / n;
                double min_x = 19.45;
                double min_y = 19.90;
                double area = 0.0;
                double peri = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        area += q([](TinyVector<double, 2> x)
                                  { return 1.0; });
                        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        peri += qp([](TinyVector<double, 2> x)
                                   { return 1.0; });
                        for (const auto &pt : q.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            qarea.evalIntegrand(xp, pt.w);
                        }
                        for (const auto &pt : qp.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            qperi.evalIntegrand(xp, pt.w);
                        }

                        for (const auto &pt : qp.nodes)
                        {
                            file << pt.x(0) << " " << pt.x(1) << " \n";
                        }
                    }
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                std::cout << "  computed area = " << area << "\n";
                double area_err = abs(area - exact_area);
                std::cout << "  area error = " << area_err << "\n";
                std::cout << "  computed perimeter = " << peri << "\n";
                double peri_err = abs(peri - 2.03955);
                std::cout << "  perimeter error = " << peri_err << "\n";
                // member function on the duration object
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                file.close();
                std::ofstream farea("airfoil_area_quadrature_rho_10.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
                std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                          << " quadrature points\n";
                std::ofstream fperi("airfoil_peri_quadrature_rho_10.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
                std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                          << " quadrature points\n";
            }
#endif
            {
                double area, peri;
                cout << "area of an airfoil using automatic subdivision " << endl;
                blitz::TinyVector<double, N> xupper;
                blitz::TinyVector<double, N> xlower;
                xlower = {19.89, 20.51};
                xupper = {19.9, 20.1};
                auto area_start = high_resolution_clock::now();
                auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, qo1);
                area = q([](const TinyVector<double, 2> x)
                         { return 1.0; });
                auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, qo1);
                peri += qp([](TinyVector<double, 2> x)
                           { return 1.0; });
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                // member function on the duration object
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                for (const auto &pt : q.nodes)
                {
                    TinyVector<double, N> xp;
                    xp[0] = pt.x[0];
                    xp[1] = pt.x[1];
                    qarea_auto.evalIntegrand(xp, pt.w);
                }
                for (const auto &pt : qp.nodes)
                {
                    TinyVector<double, N> xp;
                    xp[0] = pt.x[0];
                    xp[1] = pt.x[1];
                    qperi_auto.evalIntegrand(xp, pt.w);
                }
                double area_err = abs(area - exact_area);
                std::cout << "  area error = " << area_err << "\n";
                cout << " ----------------------- " << endl;
                std::cout << "int rule size automatic subdivision " << qarea_auto.nodes.size() << std::endl;
                std::ofstream farea("airfoil_area_quadrature_rho_10_auto.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qarea_auto, farea);
                std::cout << "  scheme.vtp file written, containing " << qarea_auto.nodes.size()
                          << " quadrature points\n";
                std::ofstream fperi("airfoil_peri_quadrature_rho_10_auto.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qperi_auto, fperi);
                std::cout << "  scheme.vtp file written, containing " << qperi_auto.nodes.size()
                          << " quadrature points\n";
            }

        } /// loop over grid size ends
    }
    if (airfoil_up)
    {
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 10;
        const int count = 1;
        double exact_area = 0.0817073;
        const int npts = 100;
        const int nbnd = npts;
        for (int q0 = 0; q0 < count; ++q0)
        {
            std::vector<TinyVector<double, N>> Xc;
            std::vector<TinyVector<double, N>> nor;
            std::vector<TinyVector<double, N - 1>> kappa;
            /// calculate boundary coordinates
            double tc = 0.12;
            TinyVector<double, npts - 1> beta;
            double beta_max = M_PI / 1.022;
            double dbeta = beta_max / (npts - 2);
            beta(0) = beta_max;
            for (int i = 1; i < npts - 1; ++i)
            {
                beta(i) = beta(i - 1) - dbeta;
            }
            TinyVector<double, nbnd> xb;
            TinyVector<double, nbnd> yb;
            constexpr bool cte = true;
            double a0 = 0.2969;
            double a1 = -0.126;
            double a2 = -0.3516;
            double a3 = 0.2843;
            double a4 = -0.1015;
            if (cte)
            {
                a4 = -0.1036;
            }
            /// upper boundary
            for (int i = 0; i < npts; ++i)
            {
                if (i == 0)
                {
                    xb(i) = 1.0;
                }
                else
                {
                    xb(i) = (1.0 - cos(beta(i - 1))) / 2.0;
                }
                double term1 = (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
                double term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
                yb(i) = 5.0 * tc * (term1 + term2);
            }

            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x;
                x(0) = xb(i);
                x(1) = yb(i);
                Xc.push_back(x);
            }
            /// get the number of boundary points
            // nbnd = Xc.size();
            cout << "nbnd " << nbnd << endl;
            cout << "Xc.size() " << Xc.size() << endl;
            /// parameters
            int ratio = 10;
            double rho = ratio * nbnd;
            cout << "rho " << rho << endl;
            /// construct the normal vector for all boundary points
            nor = constructNormal<N>(Xc);
            kappa = getCurvature<N, cte>(Xc);
            // cout << "Airfoil coordinates: " << endl;
            // for (int k = 0; k < nbnd; ++k)
            // {
            //     cout << Xc.at(k) << endl;
            // }
            /// grid size
            int n = nel[0];
            /// translate airfoil
            TinyVector<double, N> xcent;
            xcent(0) = 19.5;
            xcent(1) = 20.0;
            std::vector<TinyVector<double, N>> Xcoord;
            for (int k = 0; k < nbnd; ++k)
            {
                TinyVector<double, N> xs;
                for (int d = 0; d < N; ++d)
                {
                    xs(d) = Xc.at(k)(d) + xcent(d);
                }
                Xcoord.push_back(xs);
            }
            /// evaluate levelset and it's gradient
            Algoim::LevelSet<N> phi;
            phi.initializeLevelSet(Xcoord, nor, kappa, rho, 1.0, delta);
            TinyVector<double, N> x;
            x(0) = 19.5;
            x(1) = 20.0;
            cout << "phi " << phi(x) << endl;
            cout << "grad phi " << phi.grad(x) << endl;
            cout << "hessian phi " << phi.hessian(x) << endl;
            cout << " ------------------------------------------- " << endl;
            cout << "qo " << qo1 << endl;
            QuadratureRule<N> qarea, qperi, qarea_auto, qperi_auto;
#if 1
            /// Area and Perimeter of a 2D airfoil, computed via the cells of a Cartesian grid
            {
                const char *surface_quad_file = "quad_points_peri_auto.dat";
                ofstream file;
                file.open(surface_quad_file);
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D airfoil, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double dx = 1.0 / n;
                double dy = 0.1 / n;
                double min_x = 19.5;
                double min_y = 20.00;
                double area = 0.0;
                double peri = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        area += q([](TinyVector<double, 2> x)
                                  { return 1.0; });
                        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        peri += qp([](TinyVector<double, 2> x)
                                   { return 1.0; });
                        for (const auto &pt : q.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            qarea.evalIntegrand(xp, pt.w);
                        }
                        for (const auto &pt : qp.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            qperi.evalIntegrand(xp, pt.w);
                        }

                        for (const auto &pt : qp.nodes)
                        {
                            file << pt.x(0) << " " << pt.x(1) << " \n";
                        }
                    }
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                std::cout << "  computed area = " << area << "\n";
                double area_err = abs(area - exact_area);
                std::cout << "  area error = " << area_err << "\n";
                std::cout << "  computed perimeter = " << peri << "\n";
                double peri_err = abs(peri - 2.03955);
                std::cout << "  perimeter error = " << peri_err << "\n";
                // member function on the duration object
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                file.close();
                std::ofstream farea("airfoil_area_quadrature_rho_10.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
                std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                          << " quadrature points\n";
                std::ofstream fperi("airfoil_peri_quadrature_rho_10.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
                std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                          << " quadrature points\n";
            }
#endif
            {
                double area, peri;
                cout << "area of an airfoil using automatic subdivision " << endl;
                blitz::TinyVector<double, N> xupper;
                blitz::TinyVector<double, N> xlower;
                xlower = {19.45, 20.55};
                xupper = {19.9, 20.1};
                auto area_start = high_resolution_clock::now();
                auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, qo1);
                area = q([](const TinyVector<double, 2> x)
                         { return 1.0; });
                auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, qo1);
                peri += qp([](TinyVector<double, 2> x)
                           { return 1.0; });
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                // member function on the duration object
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                for (const auto &pt : q.nodes)
                {
                    TinyVector<double, N> xp;
                    xp[0] = pt.x[0];
                    xp[1] = pt.x[1];
                    qarea_auto.evalIntegrand(xp, pt.w);
                }
                for (const auto &pt : qp.nodes)
                {
                    TinyVector<double, N> xp;
                    xp[0] = pt.x[0];
                    xp[1] = pt.x[1];
                    qperi_auto.evalIntegrand(xp, pt.w);
                }
                double area_err = abs(area - exact_area);
                std::cout << "  area error = " << area_err << "\n";
                cout << " ----------------------- " << endl;
                std::cout << "int rule size automatic subdivision " << qarea_auto.nodes.size() << std::endl;
                std::ofstream farea("airfoil_area_quadrature_rho_10_auto.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qarea_auto, farea);
                std::cout << "  scheme.vtp file written, containing " << qarea_auto.nodes.size()
                          << " quadrature points\n";
                std::ofstream fperi("airfoil_peri_quadrature_rho_10_auto.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qperi_auto, fperi);
                std::cout << "  scheme.vtp file written, containing " << qperi_auto.nodes.size()
                          << " quadrature points\n";
            }
        }
    }
    return 0;
} // main