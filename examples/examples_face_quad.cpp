// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// on multi-component domains implicitly-defined by (one or more) multivariate Bernstein
// polynomials. Additional examples are provided on the GitHub documentation page,
// https://algoim.github.io/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "quadrature_multipoly.hpp"

using namespace algoim;

// Given a set of quadrature points and weights, output them to an VTP XML file for visualisation
// purposes, e.g., using ParaView
template <int N>
void outputQuadratureRuleAsVtpXML(const std::vector<uvector<real, N + 1>> &q, std::string fn)
{
    static_assert(N == 2 || N == 3, "outputQuadratureRuleAsVtpXML only supports 2D and 3D quadrature schemes");
    std::ofstream stream(fn);
    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "<PolyData>\n";
    stream << "<Piece NumberOfPoints=\"" << q.size() << "\" NumberOfVerts=\"" << q.size() << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    stream << "<Points>\n";
    stream << "  <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">";
    for (const auto &pt : q)
        stream << pt(0) << ' ' << pt(1) << ' ' << (N == 3 ? pt(2) : 0.0) << ' ';
    stream << "</DataArray>\n";
    stream << "</Points>\n";
    stream << "<Verts>\n";
    stream << "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
    for (size_t i = 0; i < q.size(); ++i)
        stream << i << ' ';
    stream << "</DataArray>\n";
    stream << "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
    for (size_t i = 1; i <= q.size(); ++i)
        stream << i << ' ';
    stream << "</DataArray>\n";
    stream << "</Verts>\n";
    stream << "<PointData Scalars=\"w\">\n";
    stream << "  <DataArray type=\"Float32\" Name=\"w\" NumberOfComponents=\"1\" format=\"ascii\">";
    for (const auto &pt : q)
        stream << pt(N) << ' ';
    stream << "</DataArray>\n";
    stream << "</PointData>\n";
    stream << "</Piece>\n";
    stream << "</PolyData>\n";
    stream << "</VTKFile>\n";
};


/// get volume and surface quadrature rule on cut elements
template <int N, typename F>
void GetCutInterfaceQuadScheme(const F &fphi,
                               uvector<real, N> xmin,
                               uvector<real, N> xmax,
                               const uvector<int, N> &P,
                               int q,
                               std::vector<uvector<real, N + 1>> &phase0)
{
    // Construct phi_local by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phi_local(nullptr, P);
    algoim_spark_alloc(real, phi_local);
    bernstein::bernsteinInterpolate<N>(
        [&](const uvector<real, N> &x)
        { return fphi(xmin + x * (xmax - xmin)); },
        phi_local);
    // Build quadrature hierarchy
    ImplicitPolyQuadrature<N> ipquad(phi_local);
    // Compute quadrature scheme and record the nodes & weights; phase0
    // corresponds to {phi_local < 0}
    ipquad.integrate(AlwaysGL,
                     q,
                     [&](const uvector<real, N> &x, real w)
                     {
                         if (bernstein::evalBernsteinPoly(phi_local, x) < 0)
                         {
                             phase0.push_back(add_component(x, N, w));
                         }
                     });
};

#if ALGOIM_EXAMPLES_DRIVER == 0 || ALGOIM_EXAMPLES_DRIVER == 4
int main(int argc, char *argv[])
{
#if 1
    // q-convergence study for a 2D circle
    {
        auto circle = [](const uvector<real, 2> &x)
        {
            return x(0) * x(0) + x(1) * x(1) - 4.0;
        };
        auto integrand = [](const uvector<real, 2> &x)
        {
            return 1.0;
        };
        real volume_exact = algoim::util::pi * 4.0;
        real surf_exact = 2.0 * algoim::util::pi * 2.0;
        // std::cout << "\n\nCircle q-convergence test\n";
        // std::cout << "q      area(q)         perim(q)        area error       perim error\n";
        //
        int nel = 4;
        //std::cout << "Area of a 2D circle, computed via the cells of a " << nel << " by " << nel << " Cartesian grid:\n";
        real dx = 4.2 / nel;
        const int N = 2;
        std::vector<uvector<real, N + 1>> qface_all;
        real area = 0.0;
        real peri = 0.0;
        for (int i = 0; i < nel; ++i)
            for (int j = 0; j < nel; ++j)
            {
                std::vector<uvector<real, N + 1>> surf, phase0, phase1;
                real area_c;
                real peri_c;
                uvector<real, N> xmin;
                uvector<real, N> xmax;
                xmin(0) = -2.1 + i * dx;
                xmin(1) = -2.1 + j * dx;
                xmax(0) = -2.1 + i * dx + dx;
                xmax(1) = -2.1 + j * dx + dx;
                std::vector<uvector<real, N>> qface;
                auto LSF = [&](const uvector<double, 1> &x)
                {
                    uvector<double, N> xs;
                    xs(0) = x(0);
                    xs(1) = xmax(1);
                    return circle(xs);
                };
                GetCutInterfaceQuadScheme<1>(LSF, xmin(0), xmax(0), 10, 10, qface);
                std::vector<uvector<real, N + 1>> interface_quad;
                for (int nq = 0; nq < qface.size(); ++nq)
                {
                    // std::cout << "qface.size() " << qface.size() << std::endl;
                    // std::cout << "Reference Quad rule: " << qface.at(nq) << std::endl;
                    real xq = qface.at(nq)(0);
                    uvector<real, N + 1> face_quad;
                    face_quad(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    face_quad(1) = xmax(1);
                    face_quad(2) = qface.at(nq)(1);
                    interface_quad.push_back(face_quad);
                    // std::cout << "Quad rule: " << interface_quad.at(nq) << std::endl;
                    qface_all.push_back(interface_quad.at(nq));
                }
            }
        outputQuadratureRuleAsVtpXML<N>(qface_all, "circle-interface.vtp");
    }
#endif
#if 0
    // q-convergence study for a 2D circle
    {
        auto integrand = [](const uvector<real, 1> &x)
        {
            return 1.0;
        };
        real volume_exact = algoim::util::pi * 4.0;
        real surf_exact = 2.0 * algoim::util::pi * 2.0;
        std::cout << "\n\nCircle q-convergence test\n";
        std::cout << "q      area(q)         perim(q)        area error       perim error\n";
        //
        int nel = 4;
        std::cout << "Area of a 2D circle, computed via the cells of a " << nel << " by " << nel << " Cartesian grid:\n";
        real dx = 4.2 / nel;
        const int N = 1;
        std::vector<uvector<real, N + 2>> phaseX_all, phaseY_all;
        real area = 0.0;
        real peri = 0.0;
        for (int i = 0; i < nel; ++i)
            for (int j = 0; j < nel; ++j)
            {
                std::vector<uvector<real, N + 1>> phaseX, phaseY;
                real area_c;
                real peri_c;
                uvector<real, N> xmin;
                uvector<real, N> xmax;
                xmin(0) = -2.1 + i * dx;
                xmin(1) = -2.1 + j * dx;
                xmax(0) = -2.1 + i * dx + dx;
                xmax(1) = -2.1 + j * dx + dx;
                std::cout << "xmin(0) " << xmin(0) << std::endl;
                GetFaceQuadSchemeX<1>(xmin(0), xmin(1), xmax(1), 10, 10,
                                     phaseY);
                GetFaceQuadSchemeY<1>(xmin(1), xmin(0), xmax(0), 10, 10,
                                     phaseX);
                std::vector<uvector<real, N + 2>> interface_quadY;
                for (int nq = 0; nq < phaseY.size(); ++nq)
                {

                    std::cout << "phaseY.size() " << phaseY.size() << std::endl;
                    std::cout << "Reference Quad rule: " << phaseY.at(nq) << std::endl;
                    real yq = phaseY.at(nq)(0);
                    uvector<real, N + 2> face_quad;
                    face_quad(0) = xmin(0);
                    face_quad(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    face_quad(2) = phaseY.at(nq)(1);
                    interface_quadY.push_back(face_quad);
                    std::cout << "Quad rule: " << interface_quadY.at(nq) << std::endl;
                    phaseY_all.push_back(interface_quadY.at(nq));
                }
                std::vector<uvector<real, N + 2>> interface_quadX;
                for (int nq = 0; nq < phaseX.size(); ++nq)
                {
                    std::cout << "phaseX.size() " << phaseX.size() << std::endl;
                    std::cout << "Reference Quad rule: " << phaseX.at(nq) << std::endl;
                    real xq = phaseX.at(nq)(0);
                    uvector<real, N + 2> face_quad;
                    face_quad(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    face_quad(1) = xmin(1);
                    face_quad(2) = phaseX.at(nq)(1);
                    interface_quadX.push_back(face_quad);
                    std::cout << "Quad rule: " << interface_quadX.at(nq) << std::endl;
                    phaseX_all.push_back(interface_quadX.at(nq));
                }
            }
        std::cout << "circle area: " << area << std::endl;
        std::cout << "circle perimeter: " << peri << std::endl;
        outputQuadratureRuleAsVtpXML<N + 1>(phaseX_all, "circle-phaseX.vtp");
        outputQuadratureRuleAsVtpXML<N + 1>(phaseY_all, "circle-phaseY.vtp");
    }
#endif
return 0;
}
#endif
