// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// for implicitly defined domains in hyperrectangles. The file contains a single main() routine;
// compile it as you would for any .cpp file with a main() entry point.

#include <fstream>
#include "algoim_quad.hpp"
using namespace::blitz;
template<int N>
struct Ellipsoid
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return x(0)*x(0) + 4.0*x(1)*x(1) - 1.0;
        else
            return x(0)*x(0) + 4.0*x(1)*x(1) + 9.0*x(2)*x(2) - 1.0;
    }

    template<typename T>
    blitz::TinyVector<T,N> grad(const blitz::TinyVector<T,N>& x) const
    {
        if (N == 2)
            return blitz::TinyVector<T,N>(2.0*x(0), 8.0*x(1));
        else
            return blitz::TinyVector<T,N>(2.0*x(0), 8.0*x(1), 18.0*x(2));
    }
};

template <int N>
struct Circle
{
    template <typename T>
    T operator()(const blitz::TinyVector<T, N> &x) const
    {
        return ((x(0) * x(0) + x(1) * x(1)) * (x(0) * x(0) + x(1) * x(1))) - 1.0;
    }

    template <typename T>
    blitz::TinyVector<T, N> grad(const blitz::TinyVector<T, N> &x) const
    {
        return blitz::TinyVector<T, N>(4.0 * x(0) * (x(0) * x(0) + x(1) * x(1)), 4.0 * x(1) * (x(0) * x(0) + x(1) * x(1)));
    }
};

int main(int argc, char* argv[])
{
    std::cout << "Algoim Examples - High-order quadrature algorithms for implicitly defined domains\n\n";
    std::cout << std::fixed << std::setprecision(16);

    // Area of a 2D circle, computed via the cells of a Cartesian grid
    // {
    //     int n = 16;
    //     std::cout << "Area of a 2D circle, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
    //     double dx = 2.2 / n;
    //     Circle<2> phi;
    //     double area = 0.0;
    //     for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
    //     {
    //         blitz::TinyVector<double,2> xmin = {-1.1 + i*dx, -1.1 + j*dx};
    //         blitz::TinyVector<double,2> xmax = {-1.1 + i*dx + dx, -1.1 + j*dx + dx};
    //         area += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double,2>(xmin, xmax), -1, -1, 4).sumWeights();
    //     }
    //     std::cout << "  computed area = " << area << "\n";
    //     std::cout << "    exact area = " << M_PI << "\n\n";
    // }

    
    // // "Area of a 2D circle using automatic subdivision
    // {
    //     std::cout << "Area of a 2D circle using automatic subdivision:\n";
    //     Circle<2> phi;
    //     auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double,2>(-1.1, 1.1), -1, -1, 4);
    //     double area = q([](const auto& x) { return 1.0; });
    //     std::cout << "  computed area = " << area << "\n";
    //     std::cout << "   exact area = " << M_PI << "\n";
    //     std::cout << " ========================================= " << "\n\n" ;
    // }

    // // "Area of a 2D ellipse using automatic subdivision
    // {
    //     std::cout << "Area of a 2D ellipse using automatic subdivision:\n";
    //     Ellipsoid<2> phi;
    //     auto q = Algoim::quadGen<2>(phi, Algoim::BoundingBox<double,2>(-1.1, 1.1), -1, -1, 4);
    //     double area = q([](const auto& x) { return 1.0; });
    //     std::cout << "  computed area = " << area << "\n";
    //     std::cout << "    (exact area = 1.5707963267948966)\n\n";
    // }

    // // Volume of a 3D ellipsoid using automatic subdivision
    // {
    //     std::cout << "Volume of a 3D ellipsoid using automatic subdivision:\n";
    //     Ellipsoid<3> phi;
    //     auto q = Algoim::quadGen<3>(phi, Algoim::BoundingBox<double,3>(-1.1, 1.1), -1, -1, 4);
    //     double volume = q([](const auto& x) { return 1.0; });
    //     std::cout << "  computed volume = " << volume << "\n";
    //     std::cout << "    (exact volume = 0.6981317007977318)\n\n";
    // }
    
    // Area of a 2D ellipse, computed via the cells of a Cartesian grid
    {

        int qo = 3;
        const int count = 1;
        TinyVector<double, count> Avec;
        TinyVector<double, count> Pvec;
        for (int q0 = 0; q0 < count; ++q0)
        {
            cout << "qo "  << qo << endl;
            int n = 512;
            std::cout << "Area of a 2D ellipse, computed via the cells of a " << n << " by " << n << " Cartesian grid:\n";
            double dx = 2.2 / n;
            Ellipsoid<2> phi;
            double area = 0.0;
            double peri = 0.0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    blitz::TinyVector<double, 2> xmin = {-1.1 + i * dx, -1.1 + j * dx};
                    blitz::TinyVector<double, 2> xmax = {-1.1 + i * dx + dx, -1.1 + j * dx + dx};
                    area += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), -1, -1, qo).sumWeights();
                    peri += Algoim::quadGen<2>(phi, Algoim::BoundingBox<double, 2>(xmin, xmax), 2, -1, qo).sumWeights();
                }
            double exact_area = 1.5707963267948966;
            double exact_peri = 4.844224;
            // std::cout << "  computed area = " << area << "\n";
            // std::cout << "    (exact area = 1.5707963267948966)\n\n";
            // std::cout << "  computed perimeter = " << peri << "\n";
            // std::cout << "    (exact perimeter = 4.844224110)\n\n";
            // std::cout << "area err " << abs(area - exact_area) << "\n";
            // std::cout << "peri err " << abs(peri - exact_peri) << "\n";
            Avec(q0) = abs(area - exact_area);
            Pvec(q0) = abs(peri - exact_peri);
            ++qo;
        }
        cout << "Area: " << endl;
        for (int i = 0; i < count; ++i)
        {
            cout << Avec(i) << "\n";
        }
        cout << "Perimeter: " << endl;
        for (int i = 0; i < count; ++i)
        {
            cout << Pvec(i) << "\n";
        }
    }
  

    // // Surface area of a 3D ellipsoid using automatic subdivision
    // {
    //     std::cout << "Surface area of a 3D ellipsoid using automatic subdivision:\n";
    //     Ellipsoid<3> phi;
    //     auto q = Algoim::quadGen<3>(phi, Algoim::BoundingBox<double,3>(-1.1, 1.1), 3, -1, 4);
    //     double surface_area = q.sumWeights();
    //     std::cout << "  computed surface area = " << surface_area << "\n";
    //     std::cout << "    (exact surface area = 4.4008095646649703)\n\n";
    // }

    // // Visualisation of a quadrature scheme in ParaView via XML VTP file
    // {
    //     std::cout << "Visualisation of a quadrature scheme in ParaView via XML VTP file:\n";
    //     Ellipsoid<3> phi;
    //     auto q = Algoim::quadGen<3>(phi, Algoim::BoundingBox<double,3>(-1.1, 1.1), -1, -1, 2);
    //     std::ofstream f("scheme.vtp");
    //     Algoim::outputQuadratureRuleAsVtpXML(q, f);
    //     std::cout << "  scheme.vtp file written, containing " << q.nodes.size() << " quadrature points\n";
    // }

    return 0;
}