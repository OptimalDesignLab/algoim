#include <fstream>
#include "algoim_levelset.hpp"
#include <iostream>
#include <iomanip> // std::setprecision()
bool ellipse = true;
bool circle = false;
bool circ = false;
bool airfoil = false;
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
        if (i==nbnd-1)
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
    double a4 = -0.1015; // -0.1036 for closed te
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
            kappa.push_back(1.0/rle);
        }
        else if (i==0 || i==nbnd-1)
        {
            kappa.push_back(0.0);
        }
        else
        {
            kappa.push_back(1.0/roc);
        }
        
        //kappa.push_back(0.0);
    }
    return kappa;
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
        return (((x(0) - 0.5) * (x(0) - 0.5)) + ((x(1) - 0.5) * (x(1) - 0.5)) - 0.09);
    }

    template <typename T>
    blitz::TinyVector<T, N> grad(const blitz::TinyVector<T, N> &x) const
    {
        return blitz::TinyVector<T, N>(2.0 * (x(0) - 0.5), 2.0 * (x(1) - 0.5));
    }
};

template <int N>
struct unitCircle
{
    template <typename T>
    T operator()(const blitz::TinyVector<T, N> &x) const
    {
        return (((x(0) - 0.0) * (x(0) - 0.0)) + ((x(1) - 0.0) * (x(1) - 0.0)) - 0.25);
    }

    template <typename T>
    blitz::TinyVector<T, N> grad(const blitz::TinyVector<T, N> &x) const
    {
        return blitz::TinyVector<T, N>(2.0 * (x(0) - 0.0), 2.0 * (x(1) - 0.0));
    }
};
int main(int argc, char *argv[])
{
    /// dimension
    const int N = 2;
/// use this if reading data from file
#if 0
   //const char *geometry_file = "geometry_data/circle_1.dat";
   const char *geometry_file = "geometry_data/ellipse_1.dat";
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
#endif

    double delta = 1e-10;
    if (airfoil)
    {
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 1;
        const int count = 1;
        double exact_area = 0.0817073;
        for (int q0 = 0; q0 < count; ++q0)
        {
            const char *geometry_file = "geometry_data/NACA_0012_200.dat";
            ifstream file;
            file.open(geometry_file);
            std::vector<TinyVector<double, N>> Xc;
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
                Xc.push_back(x);
            }
            file.close();
            /// get the number of boundary points
            int nbnd = Xc.size();
            cout << "nbnd " << nbnd << endl;
            /// parameters
            double rho = 10 * nbnd;
            cout << "rho " << rho << endl;
            /// construct the normal vector for all boundary points
            nor = constructNormal<N>(Xc);
            kappa = getCurvature<N>(Xc);
            // cout << "normal is: " << endl;
            // for (int k = 0; k < nbnd; ++k)
            // {
            //     cout << nor.at(k) << endl;
            // }
            /// grid size
            int n = nel[1];
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
            QuadratureRule<N> qarea, qperi, qarea_auto;
            
            /// Area of a 2D airfoil, computed via the cells of a Cartesian grid
            {
                const char *surface_quad_file = "quad_points_peri_auto.dat";
                ofstream file;
                file.open(surface_quad_file);
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D airfoil, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double dx = 1.2 / n;
                double dy = 0.2 / n;
                double min_x = 19.4;
                double min_y = 19.9;
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
                
            }
            {
                double area;
                cout << "area of an airfoil using automatic subdivision " << endl;
                blitz::TinyVector<double, N> xupper;
                blitz::TinyVector<double, N> xlower;
                xlower = {1.9, 2.4};
                xupper = {3.1, 2.6};
                auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, qo1);
                area = q([](const TinyVector<double, 2> x)
                                { return 1.0; });
                for (const auto &pt : q.nodes)
                {
                    TinyVector<double, N> xp;
                    xp[0] = pt.x[0];
                    xp[1] = pt.x[1];
                    qarea_auto.evalIntegrand(xp, pt.w);
                }
                double area_err = abs(area - exact_area);
                std::cout << "  area error = " << area_err << "\n";
                 cout << " ----------------------- " << endl;
            }
            std::cout << "int rule size automatic subdivision " << qarea_auto.nodes.size() << std::endl;
            std::ofstream farea("airfoil_rho_10_nel_16_area.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
            std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                      << " quadrature points\n";
            std::ofstream fperi("airfoil_rho_10_nel_16_peri.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
            std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                      << " quadrature points\n";
        } /// loop over grid size ends
    }
    if (circle)
    {
        const char *area = "circle_area_p2.dat";
        const char *area_err = "circle_area_err_p2.dat";
        const char *perimeter = "circle_peri_p2.dat";
        const char *perimeter_err = "circle_peri_err_p2.dat";
        ofstream file_area, file_area_err, file_peri, file_peri_err;
        file_area.open(area, ios::app);
        file_area_err.open(area_err, ios::app);
        file_peri.open(perimeter, ios::app);
        file_peri_err.open(perimeter_err, ios::app);
        file_area << std::fixed << setprecision(20) << endl;
        file_area_err << std::fixed << setprecision(20) << endl;
        file_peri << std::fixed << setprecision(20) << endl;
        file_peri_err << std::fixed << setprecision(20) << endl;
        std::vector<TinyVector<double, N>> Xc;
        std::vector<TinyVector<double, N>> nor;
        std::vector<TinyVector<double, N - 1>> kappa;
        int nel[5] = {8, 16, 32, 64, 128 /*,  32, 64, 128, 256, 512, 1024, 2048*/};
        const int count = 4;
        for (int idx = 0; idx < count; ++idx)
        {
            int n = nel[idx];
            int nbnd = n * n;
            cout << "nbnd " << nbnd << endl;
            /// parameters
            double rho = nbnd;
            /// radius
            double a = 1.0;
            /// center coords
            double xc = 0.0;
            double yc = 0.0;
            for (int k = 0; k < nbnd; ++k)
            {
                double theta = k * 2.0 * M_PI / nbnd;
                TinyVector<double, N> x, nrm;
                x(0) = a * cos(theta) + xc;
                x(1) = a * sin(theta) + yc;
                nrm(0) = 2.0 * (x(0) - xc);
                nrm(1) = 2.0 * (x(1) - yc);
                double ds = mag(nrm);
                TinyVector<double, N> ni;
                ni = nrm / ds;
                Xc.push_back(x);
                nor.push_back(ni);
                kappa.push_back(0.0);
            }
            /// evaluate levelset and it's gradient
            Algoim::LevelSet<N> phi;
            phi.initializeLevelSet(Xc, nor, kappa, rho, 1.0, delta);
            // Circle<2> phi;
            phi.xscale = 1.0;
            phi.yscale = 1.0;
            phi.min_x = 0.0;
            phi.min_y = 0.0;
            TinyVector<double, N> x;
            x(0) = 0.0;
            x(1) = 0.8;
            std::cout << std::setprecision(10) << std::endl;
            cout << "phi " << phi(x) << endl;
            cout << "grad phi " << phi.grad(x) << endl;
            // cout << "hessian phi " << phi.hessian(x) << endl;
            cout << " ------------------------------------------- " << endl;
            int qo1 = 2;
            int qo2 = 2;
            QuadratureRule<N> qa;
            /// Area of a 2D circle, computed via the cells of a Cartesian grid
            {
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D circle, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double dx = 2.2 / n;
                double dy = 2.2 / n;
                double min_x = -1.1;
                double min_y = -1.1;
                double area = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        // cout << "xmin " << xmin  << "     " << "xmax " << xmax << endl;
                        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        for (const auto &pt : q.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            qa.evalIntegrand(xp, pt.w);
                        }
                        area += q([](TinyVector<double, 2> x)
                                  { return 1.0; });
                    }
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                std::cout << "  computed area = " << area << "\n";
                std::cout << "    exact area = " << M_PI * a * a << "\n";
                double area_err = abs(area - M_PI * a * a);
                std::cout << "  area error = " << area_err << "\n";
                auto peri_start = high_resolution_clock::now();
                //  std::cout << "Peri of a 2D circle, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double peri = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        // cout << "xmin " << xmin  << "     " << "xmax " << xmax << endl;
                        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        peri += qp([](TinyVector<double, 2> x)
                                   { return 1.0; });
                    }
                std::cout << "  computed perimeter = " << peri << "\n";
                std::cout << "    (exact perimeter = " << 2 * M_PI * a << "\n";
                double peri_err = abs(peri - 2 * M_PI * a);
                std::cout << "  perimeter error = " << peri_err << "\n";
                file_area << area << " ";
                file_area_err << area_err << " ";
                file_peri << peri << " ";
                file_peri_err << peri_err << " ";
                auto peri_stop = high_resolution_clock::now();
                auto peri_duration = duration_cast<seconds>(peri_stop - peri_start);

                // To get the value of duration use the count()
                // member function on the duration object
                cout << " ---- Time taken for area ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ---- Time taken for perimeter ---- " << endl;
                cout << "      " << peri_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                std::ofstream f("element_quad_rule_ls_algoim.vtp");
                Algoim::outputQuadratureRuleAsVtpXML(qa, f);
                std::cout << "  scheme.vtp file written, containing " << qa.nodes.size()
                          << " quadrature points\n";
            }
        }
        file_area << "\n";
        file_area_err << "\n";
        file_peri << "\n";
        file_peri_err << "\n";
        file_area.close();
        file_area_err.close();
        file_peri.close();
        file_peri_err.close();

#if 0
        //"Area of a 2D circle using automatic subdivision
        {
            std::cout << "Area and Perimeter of a 2D circle using automatic subdivision:\n";
            blitz::TinyVector<double, N> xupper;
            blitz::TinyVector<double, N> xlower;
            xlower = {0.0, 0.0};
            xupper = {1.0, 1.0};
            auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, qo2);
            double area = q([](const TinyVector<double, 2> x)
                            { return 1.0; });
            std::cout << "  computed area = " << area << "\n";
            //  std::cout << "    (exact area = " << M_PI << "\n";
            auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, qo2);
            double peri = qp([](const TinyVector<double, 2> x)
                             { return 1.0; });
            std::cout << "  computed perimeter = " << peri << "\n";
            std::ofstream f("auto_element_quad_rule_ls_algoim.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(q, f);
            std::cout << "  scheme.vtp file written, containing " << q.nodes.size()
                      << " quadrature points\n";
            //std::cout << "    (exact perimeter = " << 2.0 * M_PI << "\n";
        }
#endif
#if 0
        //"Area of a 2D circle using automatic subdivision
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
            std::cout << "    (exact area = " << M_PI << "\n";
            auto qp = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, 2);
            double peri = qp([](const TinyVector<double, 2> x) { return 1.0; });
            std::cout << "  computed perimeter = " << peri << "\n";
            std::cout << "    (exact perimeter = " << 2.0 * M_PI << "\n";
        }
        {
            Circle<2> phi2;
            int n = 32;
            std::cout << "Area of a 2D circle (using exact level-set), computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
            double dx = 2.2 / n;
            double dy = 2.2 / n;
            double min_x = 0.9;
            double min_y = 0.9;
            double area = 0.0;
            double peri = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                    blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                    //cout << "xmin " << xmin  << "     " << "xmax " << xmax << endl;
                    auto q = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, 4);
                    area += q([](TinyVector<double, 2> x) { return 1.0; });
                    auto qp = Algoim::quadGen<2>(phi2, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, 4);
                    peri += qp([](TinyVector<double, 2> x) { return 1.0; });
                }
            std::cout << "  computed area = " << area << "\n";
            std::cout << "    exact area = " << M_PI << "\n";
            std::cout << "  computed perimeter = " << peri << "\n";
            std::cout << "    (exact perimeter = " << 2 * M_PI << "\n";
        }
#endif
    }

    if (ellipse)
    {
        const char *area_err = "ellipse_area_error_new_lin.dat";
        const char *perimeter_err = "ellipse_peri_error_new_lin.dat";
        ofstream file_area, file_area_err, file_peri, file_peri_err;
        file_area_err.open(area_err, ios::app);
        file_peri_err.open(perimeter_err, ios::app);
        file_area << std::fixed << setprecision(20) << endl;
        file_area_err << std::fixed << setprecision(20) << endl;
        // file_peri << std::fixed << setprecision(20) << endl;
        // file_peri_err << std::fixed << setprecision(20) << endl;
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 2;
        const int count = 5;
        TinyVector<double, count> Avec;
        TinyVector<double, count> Pvec;

        for (int q0 = 0; q0 < count; ++q0)
        {
            /// grid size
            int n = nel[q0];
            /// # boundary points
            int nbnd = 4*pow(n, qo1);
            cout << "nbnd " << nbnd << endl;
            /// parameters
            double rho = 10*nbnd;
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
            /// evaluate levelset and it's gradient
            Algoim::LevelSet<N> phi;
            phi.initializeLevelSet(Xc, nor, kappa, rho, 1.0, delta);
            TinyVector<double, N> x;
            x(0) = 0.5;
            x(1) = 0.0;
            cout << "phi " << phi(x) << endl;
            cout << "grad phi " << phi.grad(x) << endl;
            cout << "hessian phi " << phi.hessian(x) << endl;
            cout << " ------------------------------------------- " << endl;
            cout << "qo " << qo1 << endl;
            QuadratureRule<N> qarea, qperi;
            QuadratureRule<N> qarea_exact, qperi_exact;
            /// Area of a 2D ellipse, computed via the cells of a Cartesian grid
            {
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double dx = 8.2 / n;
                double dy = 2.2 / n;
                double min_x = -4.1;
                double min_y = -1.1;
                double area_c = 0.0;
                double peri_c = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        // cout << "xmin: " << xmin << " xmax: " << xmax << endl;
                        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        area_c += q([](TinyVector<double, 2> x)
                                  { return 1.0; });
                        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        peri_c += qp([](TinyVector<double, 2> x)
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
                        /// using exact lsf
                        Ellipsoid<2> phi_exact;
                        auto q_exact = Algoim::quadGen<2>(phi_exact, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        for (const auto &pt : q_exact.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            double weight = pt.w;
                            qarea_exact.evalIntegrand(xp, weight);
                        }
                        auto qp_exact = Algoim::quadGen<2>(phi_exact, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        for (const auto &pt : qp_exact.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            double weight = pt.w;
                            qperi_exact.evalIntegrand(xp, weight);
                        }
                    }
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                cout << "exact_area " << exact_area << endl;
                double area = qarea.sumWeights();
                double peri = qperi.sumWeights();
                double area_exact = qarea_exact.sumWeights();
                double peri_exact = qperi_exact.sumWeights();
                std::cout << "  computed area 1 = " << area_c << "\n";
                std::cout << "  exact area = " << exact_area << "\n";
                std::cout << "  computed area = " << area << "\n";
                double area_err = abs(area_c - exact_area);
                double exact_peri = 17.156843550313663;
                std::cout << "  area error = " << area_err << "\n";
                std::cout << "  exact lsf perimeter = " << peri_exact << "\n";
                std::cout << "  computed perimeter = " << peri << "\n";
                double peri_err = abs(peri_c - peri_exact);
                std::cout << "  perimeter error = " << peri_err << "\n";
                // member function on the duration object

                file_area_err << area_err << " ";
                file_peri_err << peri_err << " ";
                // member function on the duration object
            }
            std::ofstream farea("ellipse_quad_area.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
            std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                      << " quadrature points\n";
            std::ofstream fperi("ellipse_quad_peri.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
            std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                      << " quadrature points\n";
            // ++qo;
        } /// loop over grid size ends
        file_area_err << "\n";
        file_peri_err << "\n";
        file_area_err.close();
        file_peri_err.close();
    }

    if (circ)
    {
        const char *area_err = "circle_area_error_thesis.dat";
        const char *perimeter_err = "circle_peri_error_thesis.dat";
        ofstream file_area, file_area_err, file_peri, file_peri_err;
        file_area_err.open(area_err, ios::app);
        file_peri_err.open(perimeter_err, ios::app);
        file_area << std::fixed << setprecision(20) << endl;
        file_area_err << std::fixed << setprecision(20) << endl;
        // file_peri << std::fixed << setprecision(20) << endl;
        // file_peri_err << std::fixed << setprecision(20) << endl;
        std::cout << setprecision(20) << endl;
        int nel[9] = {8, 16, 32, 64, 128, 256, 512, 1024, 2048};
        int qo1 = 2;
        const int count = 5;
        TinyVector<double, count> Avec;
        TinyVector<double, count> Pvec;

        for (int q0 = 0; q0 < count; ++q0)
        {
            /// grid size
            int n = nel[q0];
            /// # boundary points
            int nbnd = 4*n;
            cout << "nbnd " << nbnd << endl;
            /// parameters
            double rho = 10*nbnd;
            cout << "rho " << rho << endl;
            /// major axis
            double a = 0.5;
            /// minor axis
            double b = 0.5;
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
                kappa.push_back(curv);
                //kappa.push_back(0.0);
            }
            /// evaluate levelset and it's gradient
            Algoim::LevelSet<N> phi;
            phi.initializeLevelSet(Xc, nor, kappa, rho, 1.0, delta);
            TinyVector<double, N> x;
            x(0) = 0.5;
            x(1) = 0.0;
            cout << "phi " << phi(x) << endl;
            cout << "grad phi " << phi.grad(x) << endl;
            cout << "hessian phi " << phi.hessian(x) << endl;
            cout << " ------------------------------------------- " << endl;
            cout << "qo " << qo1 << endl;
            QuadratureRule<N> qarea, qperi;
            QuadratureRule<N> qarea_exact, qperi_exact;
            /// Area of a 2D ellipse, computed via the cells of a Cartesian grid
            {
                auto area_start = high_resolution_clock::now();
                std::cout << "Area and Perimeter of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
                double dx = 1.2 / n;
                double dy = 1.2 / n;
                double min_x = -0.6;
                double min_y = -0.6;
                // double area = 0.0;
                // double peri = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                        blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                        // cout << "xmin: " << xmin << " xmax: " << xmax << endl;
                        auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        // area += q([](TinyVector<double, 2> x)
                        //           { return 1.0; });
                        auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        // peri += qp([](TinyVector<double, 2> x)
                        //            { return 1.0; });
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
                        /// using exact lsf
                        unitCircle<2> phi_exact;
                        auto q_exact = Algoim::quadGen<2>(phi_exact, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1);
                        for (const auto &pt : q_exact.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            double weight = pt.w;
                            qarea_exact.evalIntegrand(xp, weight);
                        }
                        auto qp_exact = Algoim::quadGen<2>(phi_exact, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1);
                        for (const auto &pt : qp_exact.nodes)
                        {
                            TinyVector<double, N> xp;
                            xp[0] = pt.x[0];
                            xp[1] = pt.x[1];
                            double weight = pt.w;
                            qperi_exact.evalIntegrand(xp, weight);
                        }
                    }
                auto area_stop = high_resolution_clock::now();
                auto area_duration = duration_cast<seconds>(area_stop - area_start);
                cout << " ---- Time taken  ---- " << endl;
                cout << "      " << area_duration.count() << "s " << endl;
                cout << " ----------------------- " << endl;
                cout << "exact_area " << exact_area << endl;
                double area = qarea.sumWeights();
                double peri = qperi.sumWeights();
                double area_exact = qarea_exact.sumWeights();
                double peri_exact = qperi_exact.sumWeights();
                double exact_area =  M_PI * a * a;
                std::cout << "  exact  area = " << exact_area << "\n";
                std::cout << "  computed area = " << area << "\n";
                double area_err = abs(area - exact_area);
                double exact_peri = 2.0 * M_PI * a;
                std::cout << "  area error = " << area_err << "\n";
                std::cout << "  exact  perimeter = " << exact_peri << "\n";
                std::cout << "  computed perimeter = " << peri << "\n";
                double peri_err = abs(peri - exact_peri);
                std::cout << "  perimeter error = " << peri_err << "\n";
                // member function on the duration object

                file_area_err << area_err << " ";
                file_peri_err << peri_err << " ";
                // member function on the duration object
            }
            std::ofstream farea("unitCircle_quad_area.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qarea, farea);
            std::cout << "  scheme.vtp file written, containing " << qarea.nodes.size()
                      << " quadrature points\n";
            std::ofstream fperi("unitCircle_quad_peri.vtp");
            Algoim::outputQuadratureRuleAsVtpXML(qperi, fperi);
            std::cout << "  scheme.vtp file written, containing " << qperi.nodes.size()
                      << " quadrature points\n";
            // ++qo;
        } /// loop over grid size ends
        file_area_err << "\n";
        file_peri_err << "\n";
        file_area_err.close();
        file_peri_err.close();
    }
    return 0;
} // main