
#include <fstream>
#include "algoim_quad.hpp"
using namespace blitz;
using namespace std;
using namespace Algoim;
bool ellipse = true;
bool circle = false;
template <int N>
struct Circle
{
    double operator()(const blitz::TinyVector<double, N> &x) const
    {
        return ((x(0) * x(0) + x(1) * x(1)) * (x(0) * x(0) + x(1) * x(1))) - 1.0;
    }

    blitz::TinyVector<double, N> grad(const blitz::TinyVector<double, N> &x) const
    {
        return blitz::TinyVector<double, N>(4.0 * x(0) * (x(0) * x(0) + x(1) * x(1)), 4.0 * x(1) * (x(0) * x(0) + x(1) * x(1)));
    }

    /// calculate the level-set function value (here the data type is interval)
    /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
    /// \param[out] phi - level-set function value
    Interval<N> operator()(const blitz::TinyVector<Interval<N>, N> &x) const
    {
        using std::exp;
        using std::sqrt;
        double rad = 1.0;
        /// get the centroid
        TinyVector<double, N> xc;
        xc(0) = x(0).alpha;
        xc(1) = x(1).alpha;
        double phi_xc;
        phi_xc = ((xc(0) * xc(0) + xc(1) * xc(1)) * (xc(0) * xc(0) + xc(1) * xc(1))) - 1.0;
        /// get the element domain
        TinyVector<double, N> x_max, x_min;
        x_max = xc + x(0).delta();
        x_min = xc - x(0).delta();
        TinyVector<double, N> xmax;
        xmax(0) = max(abs(x_max(0)), abs(x_min(0)));
        xmax(1) = max(abs(x_max(1)), abs(x_min(1)));
        double phi_max = ((xmax(0) * xmax(0) + xmax(1) * xmax(1)) * (xmax(0) * xmax(0) + xmax(1) * xmax(1))) - 1.0;
        /// evaluate level-set bounds
        double phi_bnd;
        phi_bnd = abs(phi_max - phi_xc);
        // cout << phi_bnd << endl;
        TinyVector<double, N> beta = grad(xc);
        double eps = phi_bnd;
        for (int dim = 0; dim < N; ++dim)
        {
            eps -= std::abs(beta(dim)) * x(dim).delta(dim);
        }
        Interval<N> phi = Interval<N>(phi_xc, beta, eps);
        return phi;
    }

    blitz::TinyVector<Interval<N>, N> grad(const blitz::TinyVector<Interval<N>, N> &x) const
    {
        using std::exp;
        using std::sqrt;
        double term1, term2;
        double rad = 1.0;
        /// get the centroid
        TinyVector<double, N> xc;
        xc(0) = x(0).alpha;
        xc(1) = x(1).alpha;
        TinyVector<double, N> delta_x, delta_y;
        delta_x = x(0).delta();
        TinyVector<double, N> phix_xc, phix_max;
        phix_xc = grad(xc);
        /// get the element domain
        TinyVector<double, N> x_max, x_min;
        x_max = xc + x(0).delta();
        x_min = xc - x(0).delta();
        double dx = x(0).delta(0);
        double dy = x(1).delta(1);
        TinyVector<double, N> xmax;
        xmax(0) = max(abs(x_max(0)), abs(x_min(0)));
        xmax(1) = max(abs(x_max(1)), abs(x_min(1)));
        phix_max = grad(xmax);
        /// evaluate level-set bounds
        // TinyVector<double, N> phix_bnd;
        // cout << phix_bnd << endl;
        TinyVector<double, N> beta_x, beta_y;
        beta_x(0) = (12.0 * xc(0) * xc(0)) + (4.0 * xc(1) * xc(1));
        beta_x(1) = 8.0 * xc(0) * xc(1);
        beta_y(0) = 8.0 * xc(0) * xc(1);
        beta_y(1) = (12.0 * xc(1) * xc(1)) + (4.0 * xc(0) * xc(0));
        TinyVector<double, N> beta_xmax, beta_ymax;
        beta_xmax(0) = (12.0 * xmax(0) * xmax(0)) + (4.0 * xmax(1) * xmax(1));
        beta_xmax(1) = 8.0 * xmax(0) * xmax(1);
        beta_ymax(0) = 8.0 * xmax(0) * xmax(1);
        beta_ymax(1) = (12.0 * xmax(1) * xmax(1)) + (4.0 * xmax(0) * xmax(0));
        term1 = beta_xmax(0) * dx;
        term2 = beta_xmax(1) * dy;
        double phix_bnd = abs(term1 + term2);
        term1 = beta_ymax(0) * dx;
        term2 = beta_ymax(1) * dy;
        double phiy_bnd = abs(term1 + term2);
        double eps_x = phix_bnd;
        double eps_y = phiy_bnd;
        // cout << phix_bnd << " , " << phiy_bnd << endl;
        // phix_bnd = abs(phix_max - phix_xc);
        for (int dim = 0; dim < N; ++dim)
        {
            eps_x -= std::abs(beta_x(dim)) * x(dim).delta(dim);
            eps_y -= std::abs(beta_y(dim)) * x(dim).delta(dim);
        }
        Interval<N> phi_x = Interval<N>(phix_xc(0), beta_x, eps_x);
        Interval<N> phi_y = Interval<N>(phix_xc(1), beta_y, eps_y);
        // Interval<N> phi_x = Interval<N>(phix_xc(0), beta_x);
        // Interval<N> phi_y = Interval<N>(phix_xc(1), beta_y);
        return blitz::TinyVector<Interval<N>, N>(phi_x, phi_y);
    }
};
template <int N>
struct Ellipse
{
    double operator()(const blitz::TinyVector<double, N> &x) const
    {
        double a = 4.0;
        double b = 1.0;
        return (((x(0) * x(0) / (a * a)) + (x(1) * x(1) / (b * b))) * ((x(0) * x(0) / (a * a)) + (x(1) * x(1) / (b * b)))) - 1.0;
    }

    blitz::TinyVector<double, N> grad(const blitz::TinyVector<double, N> &x) const
    {
        double a = 4.0;
        double b = 1.0;
        return blitz::TinyVector<double, N>(4.0 * x(0) * ((x(0) * x(0) / (a * a)) + (x(1) * x(1) / (b * b))) / (a * a),
                                            4.0 * x(1) * ((x(0) * x(0) / (a * a)) + (x(1) * x(1) / (b * b))) / (b * b));
    }

    /// calculate the level-set function value (here the data type is interval)
    /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
    /// \param[out] phi - level-set function value
    Interval<N> operator()(const blitz::TinyVector<Interval<N>, N> &x) const
    {
        using std::exp;
        using std::sqrt;
        double rad = 1.0;
        /// get the centroid
        TinyVector<double, N> xc;
        xc(0) = x(0).alpha;
        xc(1) = x(1).alpha;
        double phi_xc;
        Ellipse<N> lvst;
        phi_xc = lvst(xc);
        /// get the element domain
        TinyVector<double, N> x_max, x_min;
        x_max = xc + x(0).delta();
        x_min = xc - x(0).delta();
        TinyVector<double, N> xmax;
        xmax(0) = max(abs(x_max(0)), abs(x_min(0)));
        xmax(1) = max(abs(x_max(1)), abs(x_min(1)));
        double phi_max = lvst(xmax);
        /// evaluate level-set bounds
        double phi_bnd;
        phi_bnd = abs(phi_max - phi_xc);
        // cout << phi_bnd << endl;
        TinyVector<double, N> beta = grad(xc);
        double eps = phi_bnd;
        for (int dim = 0; dim < N; ++dim)
        {
            eps -= std::abs(beta(dim)) * x(dim).delta(dim);
        }
        Interval<N> phi = Interval<N>(phi_xc, beta, eps);
        return phi;
    }

    blitz::TinyVector<Interval<N>, N> grad(const blitz::TinyVector<Interval<N>, N> &x) const
    {
        using std::exp;
        using std::sqrt;
        double term1, term2;
        double rad = 1.0;
        double a = 4.0;
        double b = 1.0;
        /// get the centroid
        TinyVector<double, N> xc;
        xc(0) = x(0).alpha;
        xc(1) = x(1).alpha;
        TinyVector<double, N> delta_x, delta_y;
        delta_x = x(0).delta();
        TinyVector<double, N> phix_xc, phix_max;
        phix_xc = grad(xc);
        /// get the element domain
        TinyVector<double, N> x_max, x_min;
        x_max = xc + x(0).delta();
        x_min = xc - x(0).delta();
        double dx = x(0).delta(0);
        double dy = x(1).delta(1);
        TinyVector<double, N> xmax;
        xmax(0) = max(abs(x_max(0)), abs(x_min(0)));
        xmax(1) = max(abs(x_max(1)), abs(x_min(1)));
        phix_max = grad(xmax);
        /// evaluate level-set bounds
        // TinyVector<double, N> phix_bnd;
        // cout << phix_bnd << endl;
        TinyVector<double, N> beta_x, beta_y;
        double fac1 = a * a * a * a;
        double fac2 = b * b * b * b;
        beta_x(0) = (8 * xc(0) * xc(0) / fac1) +
                    (4.0 * ((xc(0) * xc(0) / (a * a)) + (xc(1) * xc(1) / (b * b))) / (a * a));
        beta_x(1) = 8.0 * xc(0) * xc(1) / (a * a * b * b);
        beta_y(0) = 8.0 * xc(0) * xc(1) / (a * a * b * b);
        beta_y(1) = (8 * xc(1) * xc(1) / fac2) +
                    (4.0 * ((xc(0) * xc(0) / (a * a)) + (xc(1) * xc(1) / (b * b))) / (b * b));
        TinyVector<double, N> beta_xmax, beta_ymax;
        beta_xmax(0) = (8 * xmax(0) * xmax(0) / fac1) +
                       (4.0 * ((xmax(0) * xmax(0) / (a * a)) + (xmax(1) * xmax(1) / (b * b))) / (a * a));
        beta_xmax(1) = 8.0 * xmax(0) * xmax(1) / (a * a * b * b);
        beta_ymax(0) = 8.0 * xmax(0) * xmax(1) / (a * a * b * b);
        beta_ymax(1) = (8 * xmax(1) * xmax(1) / fac2) +
                       (4.0 * ((xmax(0) * xmax(0) / (a * a)) + (xmax(1) * xmax(1) / (b * b))) / (b * b));
        term1 = beta_xmax(0) * dx;
        term2 = beta_xmax(1) * dy;
        double phix_bnd = abs(term1 + term2);
        term1 = beta_ymax(0) * dx;
        term2 = beta_ymax(1) * dy;
        double phiy_bnd = abs(term1 + term2);
        double eps_x = phix_bnd;
        double eps_y = phiy_bnd;
        // cout << phix_bnd << " , " << phiy_bnd << endl;
        // phix_bnd = abs(phix_max - phix_xc);
        for (int dim = 0; dim < N; ++dim)
        {
            eps_x -= std::abs(beta_x(dim)) * x(dim).delta(dim);
            eps_y -= std::abs(beta_y(dim)) * x(dim).delta(dim);
        }
        Interval<N> phi_x = Interval<N>(phix_xc(0), beta_x, eps_x);
        Interval<N> phi_y = Interval<N>(phix_xc(1), beta_y, eps_y);
        // Interval<N> phi_x = Interval<N>(phix_xc(0), beta_x);
        // Interval<N> phi_y = Interval<N>(phix_xc(1), beta_y);
        return blitz::TinyVector<Interval<N>, N>(phi_x, phi_y);
    }
};

/// This functions evaluates the level-set bounds along given x = xslice
/// \param[in] xslice - slice in x direction
/// \param[in] ymax - maximum `y` coordinate
/// \param[in] ymin - minimum `y` coordinate
/// \param[in] nel - number of elements in `y` direction
template <int N>
void testLevelSetBounds(double xslice, double ymax, double ymin, const int nel)
{
    /// create level-set object
    Ellipse<N> phi;
    TinyVector<Interval<N>, N> phic;
    double ds = (ymax - ymin) / nel;
    TinyVector<double, N> del;
    del = {0.5 * ds, 0.5 * ds};
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
/// \param[in] xslice - slice in x direction
/// \param[in] ymax - maximum `y` coordinate
/// \param[in] ymin - minimum `y` coordinate
/// \param[in] nel - number of elements in `y` direction
template <int N>
void testLevelSetGradBounds(double xslice, double ymax, double ymin, const int nel)
{
    /// create level-set object
    Ellipse<N> phi;
    TinyVector<Interval<N>, N> phic;
    double ds = (ymax - ymin) / nel;
    TinyVector<double, N> del;
    del = {0.5 * ds, 0.5 * ds};
    cout << " ------------------------------------------- " << endl;
    cout << "                  Gradient bounds      " << endl;
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

int main(int argc, char *argv[])
{
    std::cout << "Algoim Examples - High-order quadrature algorithms for implicitly defined domains\n\n";
    std::cout << std::fixed << std::setprecision(16);
    /// get the bounds
    const int nel = 20;
    double ymin = -2.0;
    double ymax = 2.0;
    double xslice = 3.0;
    cout << setprecision(12) << endl;
    const int N = 2;
    int qo1 = 10;
    int qo2 = 10;
    // testLevelSetBounds<N>(xslice, ymax, ymin, nel);
    // testLevelSetGradBounds<N>(xslice, ymax, ymin, nel);
    if (ellipse)
    {
        Ellipse<2> phi;
        // Area of a 2D ellipse, computed via the cells of a Cartesian grid
        {
            int n = 64;
            std::cout << "Area and Perimeter of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
            double dx = 8.2/ n;
            double dy = 2.2 / n;
            double min_x = -4.1;
            double min_y = -1.1;
            double area = 0.0;
            double peri = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    blitz::TinyVector<double, 2> xmin = {min_x + i * dx, min_y + j * dy};
                    blitz::TinyVector<double, 2> xmax = {min_x + i * dx + dx, min_y + j * dy + dy};
                    area += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo1).sumWeights();
                    peri += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo1).sumWeights();
                }
            std::cout << "  computed area = " << area << "\n";
            std::cout << "    exact area = " << 4.0 * M_PI << "\n\n";
            std::cout << "  computed peri = " << peri << "\n";
            std::cout << "    exact peri = " << 17.15684 << "\n\n";
        }

        // "Area of a 2D ellipse using automatic subdivision
        {
            std::cout << "Area and Perimeter of a 2D ellipse using automatic subdivision:\n";
            blitz::TinyVector<double, N> xupper;
            blitz::TinyVector<double, N> xlower;
            xlower = {-4.1, -1.1};
            xupper = {4.1, 1.1};
            auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), -1, -1, qo2);
            double area = q([](const auto &x) { return 1.0; });
            auto qp = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xlower, xupper), 2, -1, qo2);
            double peri = qp([](const auto &x) { return 1.0; });
            std::cout << "  computed area = " << area << "\n";
            std::cout << "   exact area = " << 4.0 * M_PI << "\n\n";
            std::cout << "  computed peri = " << peri << "\n";
            std::cout << "    exact peri = " << 17.15684 << "\n\n";
            std::cout << " ========================================= "
                      << "\n\n";
        }
    }
    if (circle)
    {
        // Area of a 2D circle, computed via the cells of a Cartesian grid
        {
            int n = 64;
            std::cout << "Area of a 2D circle, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
            double dx = 2.2 / n;
            Circle<2> phi;
            double area = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    blitz::TinyVector<double, 2> xmin = {-1.1 + i * dx, -1.1 + j * dx};
                    blitz::TinyVector<double, 2> xmax = {-1.1 + i * dx + dx, -1.1 + j * dx + dx};
                    area += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, 4).sumWeights();
                }
            std::cout << "  computed area = " << area << "\n";
            std::cout << "    exact area = " << M_PI << "\n\n";
        }

        // Perimeter of a 2D circle, computed via the cells of a Cartesian grid
        {
            int n = 64;
            std::cout << "Perimeter of a 2D circle, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
            double dx = 2.2 / n;
            Circle<2> phi;
            double peri = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    blitz::TinyVector<double, 2> xmin = {-1.1 + i * dx, -1.1 + j * dx};
                    blitz::TinyVector<double, 2> xmax = {-1.1 + i * dx + dx, -1.1 + j * dx + dx};
                    peri += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, 4).sumWeights();
                }
            std::cout << "  computed peri = " << peri << "\n";
            std::cout << "    exact peri = " << 2.0 * M_PI << "\n\n";
        }

        // "Area of a 2D circle using automatic subdivision
        {
            std::cout << "Area of a 2D circle using automatic subdivision:\n";
            Circle<2> phi;
            auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(-1.1, 1.1), -1, -1, 4);
            double area = q([](const auto &x) { return 1.0; });
            std::cout << "  computed area = " << area << "\n";
            std::cout << "   exact area = " << M_PI << "\n\n";
        }

        // "Perimeter of a 2D circle using automatic subdivision
        {
            std::cout << "Perimeter of a 2D circle using automatic subdivision:\n";
            Circle<2> phi;
            auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(-1.1, 1.1), 2, -1, 4);
            double peri = q([](const auto &x) { return 1.0; });
            std::cout << "  computed peri = " << peri << "\n";
            std::cout << "    exact peri = " << 2.0 * M_PI << "\n\n";
            std::cout << " ========================================= "
                      << "\n\n";
        }
    }

    return 0;
}
