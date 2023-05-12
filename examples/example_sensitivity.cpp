// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// on multi-component domains implicitly-defined by (one or more) multivariate Bernstein
// polynomials. Additional examples are provided on the GitHub documentation page,
// https://algoim.github.io/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "quadrature_multipoly.hpp"
#include "duals/dual"
using namespace algoim;
using duals::duald;
using namespace duals::literals;

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
template <typename xdouble, int N, typename F>
void GetCutInterfaceQuadScheme(const F &fphi,
                               uvector<real, N> xmin,
                               uvector<real, N> xmax,
                               const uvector<int, N> &P,
                               int q,
                               std::vector<uvector<xdouble, N + 1>> &phase0)
{
    // Construct phi_local by mapping [0,1] onto bounding box [xmin,xmax]
    std::cout << "================================================================ " << std::endl;
    std::cout << "Inside GetCutInterfaceQuadScheme ( ) " << std::endl;
    xarray<xdouble, N> phi_local(nullptr, P);
    // std::cout << "created phi_local " << std::endl;
    algoim_spark_alloc(xdouble, phi_local);
    // std::cout << "allocated phi_local " << std::endl;
    // std::cout << "before interpolation " << std::endl;
    bernstein::bernsteinInterpolate<xdouble, N>(
        [&](const uvector<real, N> &x)
        { return fphi(xmin + x * (xmax - xmin)); },
        phi_local);

    for (int k = 0; k < P(0); ++k)
    {
        std::cout << "phi_local[" << k << "]: " << phi_local[k] << std::endl;
    }
    //std::cout << "interpolation is done !! " << std::endl;
    // Build quadrature hierarchy
    ImplicitPolyQuadrature<xdouble, N> ipquad(phi_local);
    // Compute quadrature scheme and record the nodes & weights; phase0
    // corresponds to {phi_local < 0}
    ipquad.integrate(AutoMixed,
                     q,
                     [&](const uvector<xdouble, N> &x, xdouble w)
                     {
                         if (bernstein::evalBernsteinPoly(phi_local, x) < 0)
                         {
                            //  for (int dim = 0; dim < N; ++dim)
                            //  {
                            //      std::cout << "x " << x(dim) << std::endl;
                            //      std::cout << "w " << w << std::endl;
                            //  }
                             phase0.push_back(add_component(x, N, w));
                         }
                     });
    std::cout << "quad size: " << phase0.size() << std::endl;
};
/// get volume and surface quadrature rule on cut elements
template <typename xdouble, int N, typename F>
void GetImmersedBoundaryQuadScheme(const F &fphi,
                                   uvector<real, N> xmin,
                                   uvector<real, N> xmax,
                                   const uvector<int, N> &P,
                                   int q,
                                   std::vector<uvector<xdouble, N + 1>> &surf)
{
    // Construct phi_local by mapping [0,1] onto bounding box [xmin,xmax]
    //
    std::cout << "Inside GetCutInterfaceQuadScheme ( ) " << std::endl;
    xarray<xdouble, N> phi_local(nullptr, P);
    std::cout << "created phi_local " << std::endl;
    algoim_spark_alloc(xdouble, phi_local);
    std::cout << "allocated phi_local " << std::endl;
    std::cout << "before interpolation " << std::endl;
    bernstein::bernsteinInterpolate<xdouble, N>(
        [&](const uvector<real, N> &x)
        { return fphi(xmin + x * (xmax - xmin)); },
        phi_local);

    // for (int k = 0; k < P(0); ++k)
    // {
    //     std::cout << "phi_local[" << k << "]: " << phi_local[k] << std::endl;
    // }
    // std::cout << "interpolation is done !! " << std::endl;
    // Build quadrature hierarchy
    ImplicitPolyQuadrature<xdouble, N> ipquad(phi_local);
    // Compute quadrature scheme and record the nodes & weights; phase0
    // corresponds to {phi_local < 0}
    ipquad.integrate_surf(AutoMixed, q,
                          [&](const uvector<xdouble, N> &x, xdouble w, const uvector<xdouble, N> &wn)
                          {
                              surf.push_back(add_component(x, N, w));
                          });
    std::cout << "surface quad size: " << surf.size() << std::endl;
    std::cout << "================================================================ " << std::endl;
};
/// get volume and surface quadrature rule on cut elements
template <typename xdouble, int N>
void computeDesignSensitivity(xdouble r_a, xdouble *w_a, xdouble *x_a)
{
    auto circle = [&](const uvector<xdouble, 2> &x)
    {
        return x(0) * x(0) + x(1) * x(1) - r_a * r_a;
    };
    std::cout << "r = " << r_a << std::endl;
    /// bounding box/line
    uvector<real, N> xmin_face;
    uvector<real, N> xmax_face;
    xmin_face(0) = 1.05;
    xmin_face(1) = 0.0;
    xmax_face(0) = 2.1;
    xmax_face(1) = 0.0;
    /// bounding box (immersed surface)
    uvector<real, N> xmin;
    uvector<real, N> xmax;
    xmin(0) = -2.1;
    xmin(1) = -2.1;
    xmax(0) = 2.1;
    xmax(1) = 2.1;
    /// general level-set function restricted to 1-D
    auto levelSetFace = [&](const uvector<real, 1> &x)
    {
        uvector<xdouble, N> xs;
        xs(0) = x(0);
        xs(1) = xmax_face(1);
        xdouble y = circle(xs);
        return y;
    };
    /// general level-set function restricted to 1-D
    auto LSF = [&](const uvector<real, 2> &x)
    {
        uvector<xdouble, N> xs;
        xs(0) = x(0);
        xs(1) = x(1);
        xdouble y = circle(xs);
        return y;
    };
    uvector<real, N - 1> xs;
    xs(0) = 2.0;
    std::cout << "levelSetFace(x) " << levelSetFace(xs) << std::endl;
    uvector<real, N> xsurf;
    xsurf(0) = 2.0;
    xsurf(1) = 0.0;
    std::cout << "LSF(xsurf) " << LSF(xsurf) << std::endl;
    /// polynominal degree and quadrature order
    int P = 3;
    int q = 3;
    /// get the interface quadrature rule
    std::vector<uvector<xdouble, N>> qface;
    std::vector<uvector<xdouble, N + 1>> surf;
    // GetCutInterfaceQuadScheme<xdouble, 1>(levelSetFace, xmin_face(0), xmax_face(0), P, q, qface);
    GetImmersedBoundaryQuadScheme<xdouble, 2>(LSF,
                                              xmin,
                                              xmax,
                                              P,
                                              q,
                                              surf);
#if 0
    for (int nq = 0; nq < qface.size(); ++nq)
    {
        // std::cout << "nq " << nq << std::endl;
        // std::cout << qface.at(nq)(0) << std::endl;
        // x_a.at(nq) = qface.at(nq)(0);
        // w_a.at(nq) = qface.at(nq)(1);
        x_a[nq] = qface.at(nq)(0);
        w_a[nq] = qface.at(nq)(1);
        // std::cout << "quad x: " << x_a.at(nq) << std::endl;
    }
#endif
#if 1
    for (int nq = 0; nq < surf.size(); ++nq)
    {
        x_a[nq] = surf.at(nq)(0);
        w_a[nq] = surf.at(nq)(2);
    }
#endif
    // calculate the jacobian w.r.t state vaiables
/// output the quadrature rule
#if 0
    std::vector<uvector<real, N + 1>> qface_all;
    std::vector<uvector<real, N + 1>> interface_quad;
    for (int nq = 0; nq < qface.size(); ++nq)
    {
        real xq = qface.at(nq)(0);
        uvector<real, N + 1> face_quad;
        face_quad(0) = xmin(0) + xq * (xmax(0) - xmin(0));
        face_quad(1) = xmax(1);
        face_quad(2) = qface.at(nq)(1);
        interface_quad.push_back(face_quad);
        qface_all.push_back(interface_quad.at(nq));
    }
    std::cout << "************************************************************************ " << std::endl;
    outputQuadratureRuleAsVtpXML<N>(qface_all, "circle-single-interface.vtp");
#endif
#if 1
    std::vector<uvector<real, N + 1>> qsurf_all;
    std::vector<uvector<real, N + 1>> surface_quad;
    if constexpr (std::is_same_v<double, xdouble>)
    {
        for (int nq = 0; nq < surf.size(); ++nq)
        {
            real xq = surf.at(nq)(0);
            real yq = surf.at(nq)(1);
            uvector<real, N + 1> surf_quad;
            surf_quad(0) = xmin(0) + xq * (xmax(0) - xmin(0));
            surf_quad(1) = xmin(1) + yq * (xmax(1) - xmin(1));
            surf_quad(2) = surf.at(nq)(2);
            surface_quad.push_back(surf_quad);
        }
    }
    else
    {
        for (int nq = 0; nq < surf.size(); ++nq)
        {
            real xq = surf.at(nq)(0).rpart();
            real yq = surf.at(nq)(1).rpart();
            uvector<real, N + 1> surf_quad;
            surf_quad(0) = xmin(0) + xq * (xmax(0) - xmin(0));
            surf_quad(1) = xmin(1) + yq * (xmax(1) - xmin(1));
            surf_quad(2) = surf.at(nq)(2).rpart();
            surface_quad.push_back(surf_quad);
        }
    }
    std::cout << "************************************************************************ " << std::endl;
    outputQuadratureRuleAsVtpXML<N>(surface_quad, "circle-boundary.vtp");
#endif
}

#if ALGOIM_EXAMPLES_DRIVER == 0 || ALGOIM_EXAMPLES_DRIVER == 4
int main(int argc, char *argv[])
{
    // real dwdr = 0.0;
    std::cout << std::setprecision(12);
    const int N = 2;
    int qs = 12;
#if 1
    duald r = 2.0 + 1_e;
    std::vector<duald> w_a(qs);
    std::vector<duald> x_a(qs);
    computeDesignSensitivity<duald, N>(r, w_a.data(), x_a.data());
    //     std::cout << "+++++++++++++++++++++++++++++++++++++ " << std::endl;
    // std::cout << "dxdr: " << std::endl;
    // for (int k = 0; k < qs; ++k)
    // {
    //     std::cout << x_a.at(k).dpart() << " ";
    // }
    // std::cout << "\n"
    //           << "+++++++++++++++++++++++++++++++++++++ " << std::endl;
#endif
#if 0
    real r = 2.0 ;
    real w;
    std::vector<real> w_a(qs);
    std::vector<real> x_a(qs);
    computeDesignSensitivity<real, N>(r, w_a.data(), x_a.data());
#endif

#if 1
    real rr = 2.00;
    real wp, wm;
    real delta = 1e-03;
    std::vector<real> w_ap(qs);
    std::vector<real> w_am(qs);
    std::vector<real> x_ap(qs);
    std::vector<real> x_am(qs);
    computeDesignSensitivity<real, N>(rr + delta, w_ap.data(), x_ap.data());
    computeDesignSensitivity<real, N>(rr - delta, w_am.data(), x_am.data());
    std::vector<real> dxdr(qs);
    std::vector<real> dydr(qs);
    std::vector<real> dwdr(qs);
    for (int k = 0; k < qs; ++k)
    {
        dwdr.at(k) = (w_ap.at(k) - w_am.at(k)) / (2.0 * delta);
        dxdr.at(k) = (x_ap.at(k) - x_am.at(k)) / (2.0 * delta);
    }
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    std::cout << "dxdr: " << std::endl;
    for (int k = 0; k < qs; ++k)
    {
        std::cout << x_a.at(k).dpart() << " ";
    }
    std::cout << "\n"
              << "dxdr-FD: " << std::endl;
    for (int k = 0; k < qs; ++k)
    {
        std::cout << dxdr.at(k) << " ";
    }
    std::cout << "\n"
              << "dwdr: " << std::endl;
    for (int k = 0; k < qs; ++k)
    {
        std::cout << w_a.at(k).dpart() << " ";
    }
    std::cout << "\n"
              << "dwdr-FD: " << std::endl;
    for (int k = 0; k < qs; ++k)
    {
        std::cout << dwdr.at(k) << " ";
    }
    std::cout << "\n"
              << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    std::cout << "abs(dxdr-duals - dxdr-FD) : " << std::endl;
    for (int k = 0; k < qs; ++k)
    {
        std::cout << abs(dxdr.at(k) - x_a.at(k).dpart()) << " ";
    }
    std::cout << "\n";
    std::cout << "abs(dwdr-duals - dwdr-FD) : " << std::endl;
    for (int k = 0; k < qs; ++k)
    {
        std::cout << abs(dwdr.at(k) - w_a.at(k).dpart()) << " ";
    }
    std::cout << "\n";
    std::cout << "************************************************************************************************************  " << std::endl;
#endif

#if 0
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
        // std::cout << "Area of a 2D circle, computed via the cells of a " << nel << " by " << nel << " Cartesian grid:\n";
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
                auto LSF = [&](const uvector<real, 1> &x)
                {
                    uvector<real, N> xs;
                    xs(0) = x(0);
                    xs(1) = xmax(1);
                    return circle(xs);
                };
                GetCutInterfaceQuadScheme<adouble, 1>(LSF, xmin(0), xmax(0), 10, 10, qface);
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
        std::cout << "************************************************************************ " << std::endl;
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
