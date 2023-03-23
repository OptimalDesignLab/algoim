// Examples to demonstrate Algoim's methods for computing high-order accurate quadrature schemes
// on multi-component domains implicitly-defined by (one or more) multivariate Bernstein
// polynomials. Additional examples are provided on the GitHub documentation page,
// https://algoim.github.io/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "quadrature_multipoly.hpp"

using namespace algoim;
/// this function constructs the normal vectors for given boundary points of level-set geometry
template <int N>
std::vector<uvector<real, N>> constructNormal(std::vector<uvector<real, N>> Xc)
{
    std::vector<uvector<real, N>> nor;
    int nbnd = Xc.size();
    uvector<real, N> nsurf;
    for (int i = 0; i < nbnd; i++)
    {
        std::vector<real> dX;
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
        real nx = dX.at(1);
        real ny = -dX.at(0);
        real ds = 1.0 / sqrt((nx * nx) + (ny * ny));
        nsurf(0) = nx * ds;
        nsurf(1) = ny * ds;
        // if (i == 0 || i == nbnd - 1)
        // {
        //     nsurf(0) = 1.0;
        //     nsurf(1) = 0;
        // }
        // if (i == 0)
        // {
        //     nsurf(0) = 1.0;
        //     nsurf(1) = 0;
        // }
        nor.push_back(nsurf);
    }
    return nor;
}
/// this function constructs the curvature at given boundary points of level-set geometry
template <int N, bool cte = false>
std::vector<uvector<real, N - 1>> getCurvature(std::vector<uvector<real, N>> Xc)
{
    std::vector<uvector<real, N - 1>> kappa;
    int nbnd = Xc.size();
    real a0 = 0.2969;
    real a1 = -0.126;
    real a2 = -0.3516;
    real a3 = 0.2843;
    real a4 = -0.1015;
    if (cte)
    {
        a4 = -0.1036;
    }
    real tc = 0.12;
    real theta = 0.0;
    // thickness
    for (int i = 0; i < nbnd; i++)
    {
        real xc = Xc.at(i)(0);
        real yt = tc * ((a0 * sqrt(xc)) + (a1 * xc) + (a2 * (pow(xc, 2))) + (a3 * (pow(xc, 3))) + (a4 * (pow(xc, 4)))) / (0.2);
        real dytdx = tc * ((0.5 * a0 / sqrt(xc)) + a1 + (2.0 * a2 * xc) + (3.0 * a3 * pow(xc, 2)) + (4.0 * a4 * pow(xc, 3))) / (0.2);
        real d2ytdx = tc * ((-0.25 * a0 / pow(xc, 1.5)) + (2.0 * a2) + (6.0 * a3 * xc) + (12.0 * a4 * pow(xc, 2))) / 0.2;
        real ysu = yt * cos(theta);
        real dydx = dytdx * cos(theta);
        real d2ydx = d2ytdx * cos(theta);
        real roc = (pow((1 + (dydx * dydx)), 1.5)) / abs(d2ydx);
        if (xc == 0.0)
        {
            real rle = 0.5 * pow(a0 * tc / 0.20, 2);
            kappa.push_back(1.0 / rle);
        }
        else if (i == 0 || i == 1 || i == nbnd - 1 || i == nbnd - 2)
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

// Driver method which takes a functor phi defining a single polynomial in the reference
// rectangle [xmin, xmax]^N, of Bernstein degree P, builds a quadrature scheme with the
// given q, and outputs it for visualisation in a set of VTP XML files
template <int N, typename F>
void outputQuadScheme(const F &fphi, real xmin, real xmax, const uvector<int, N> &P, int q, std::string qfile)
{
    // Construct phi by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phi(nullptr, P);
    algoim_spark_alloc(real, phi);
    bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
                                       { return fphi(xmin + x * (xmax - xmin)); },
                                       phi);

    // Build quadrature hierarchy
    ImplicitPolyQuadrature<N> ipquad(phi);
    // ipquad.type = OuterAggregate;
    // Compute quadrature scheme and record the nodes & weights; phase0 corresponds to
    // {phi < 0}, phase1 corresponds to {phi > 0}, and surf corresponds to {phi == 0}.
    std::vector<uvector<real, N + 1>> phase0, phase1, surf;
    ipquad.integrate(AutoMixed, q, [&](const uvector<real, N> &x, real w)
                     {
        if (bernstein::evalBernsteinPoly(phi, x) < 0)
            phase0.push_back(add_component(x, N, w));
        else
            phase1.push_back(add_component(x, N, w)); });
    ipquad.integrate_surf(AutoMixed, q, [&](const uvector<real, N> &x, real w, const uvector<real, N> &wn)
                          { surf.push_back(add_component(x, N, w)); });
    // ipquad.integrate_surf(AlwaysGL, q, [&](const uvector<real, N> &x, real w, const uvector<real, N> &wn)
    //                       { surf.push_back(add_component(x, N, w)); });
    std::cout << "Surface quadrature: " << surf.size() << std::endl;
    for (int i = 0; i < surf.size(); ++i)
    {
        // std::cout << surf.at(i) << std::endl;
        real xq = surf.at(i)(0);
        real yq = surf.at(i)(1);
        real dy = 0.0;
        // surf.at(i)(0) = xmin + xq * (xmax - xmin);
        // surf.at(i)(1) = xmin + yq * (xmax - xmin) + dy;
    }
    // output to file
    outputQuadratureRuleAsVtpXML<N>(phase0, qfile + "-phase0.vtp");
    outputQuadratureRuleAsVtpXML<N>(phase1, qfile + "-phase1.vtp");
    outputQuadratureRuleAsVtpXML<N>(surf, qfile + "-surf.vtp");
}
/// perform h refinement study
template <int N, typename F>
void GetQuadScheme(const F &fphi, uvector<real, N> xmin, uvector<real, N> xmax, const uvector<int, N> &P, int q, std::vector<uvector<real, N + 1>> &surf,
                   std::vector<uvector<real, N + 1>> &phase0, std::vector<uvector<real, N + 1>> &phase1)
{
    // Construct phi by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phi(nullptr, P);
    algoim_spark_alloc(real, phi);
    bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
                                       { return fphi(xmin + x * (xmax - xmin)); },
                                       phi);
    // Build quadrature hierarchy
    ImplicitPolyQuadrature<2> ipquad(phi);
    // ipquad.type = OuterAggregate;
    // Compute quadrature scheme and record the nodes & weights; phase0 corresponds to
    // {phi < 0}, phase1 corresponds to {phi > 0}, and surf corresponds to {phi == 0}.
    ipquad.integrate(AlwaysGL, q, [&](const uvector<real, N> &x, real w)
                     {
        if (bernstein::evalBernsteinPoly(phi, x) < 0)
            phase0.push_back(add_component(x, N, w));
        else
            phase1.push_back(add_component(x, N, w)); });
    ipquad.integrate_surf(AlwaysGL, q, [&](const uvector<real, N> &x, real w, const uvector<real, N> &wn)
                          { surf.push_back(add_component(x, N, w)); });
}

// Driver method which takes a functor phi defining a single polynomial in the reference
// rectangle [xmin, xmax]^N, of Bernstein degree P, along with an integrand function,
// and performances a q-refinement convergence study, comparing the computed integral
// with the given exact answers, for 1 <= q <= qMax.
template <int N, typename Phi, typename F>
void qConv(const Phi &phi, real xmin, real xmax, uvector<int, N> P, const F &integrand, int qMax, real volume_exact, real surf_exact)
{
    // Construct Bernstein polynomial by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phipoly(nullptr, P);
    algoim_spark_alloc(real, phipoly);
    bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
                                       { return phi(xmin + x * (xmax - xmin)); },
                                       phipoly);

    // Build quadrature hierarchy
    ImplicitPolyQuadrature<N> ipquad(phipoly);

    // Functional to evaluate volume and surface integrals of given integrand
    real volume, surf;
    auto compute = [&](int q)
    {
        volume = 0.0;
        surf = 0.0;
        // compute volume integral over {phi < 0} using AutoMixed strategy
        ipquad.integrate(AutoMixed, q, [&](const uvector<real, N> &x, real w)
                         {
            if (bernstein::evalBernsteinPoly(phipoly, x) < 0)
                volume += w * integrand(xmin + x * (xmax - xmin)); });
        // compute surface integral over {phi == 0} using AutoMixed strategy
        ipquad.integrate_surf(AutoMixed, q, [&](const uvector<real, N> &x, real w, const uvector<real, N> &wn)
                              { surf += w * integrand(xmin + x * (xmax - xmin)); });
        // scale appropriately
        volume *= pow(xmax - xmin, N);
        surf *= pow(xmax - xmin, N - 1);
    };

    // Compute results for all q and output in a convergence table
    for (int q = 1; q <= qMax; ++q)
    {
        compute(q);
        std::cout << q << ' ' << volume << ' ' << surf << ' ' << std::abs(volume - volume_exact) / volume_exact << ' ' << std::abs(surf - surf_exact) / surf_exact << std::endl;
    }
}
template <int N, typename Phi, typename F>
void qAreaPeri(const Phi &phi, uvector<real, N> xmin, uvector<real, N> xmax, uvector<int, N> P, const F &integrand, int qMax, real volume_exact, real surf_exact, real &volume, real &surf)
{
    // Construct Bernstein polynomial by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phipoly(nullptr, P);
    algoim_spark_alloc(real, phipoly);
    bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
                                       { return phi(xmin + x * (xmax - xmin)); },
                                       phipoly);

    // Build quadrature hierarchy
    ImplicitPolyQuadrature<N> ipquad(phipoly);

    // Functional to evaluate volume and surface integrals of given integrand
    auto compute = [&](int q)
    {
        volume = 0.0;
        surf = 0.0;
        // compute volume integral over {phi < 0} using AutoMixed strategy
        ipquad.integrate(AutoMixed, q, [&](const uvector<real, N> &x, real w)
                         {
            if (bernstein::evalBernsteinPoly(phipoly, x) < 0)
                volume += w * integrand(xmin + x * (xmax - xmin)); });
        // compute surface integral over {phi == 0} using AutoMixed strategy
        ipquad.integrate_surf(AutoMixed, q, [&](const uvector<real, N> &x, real w, const uvector<real, N> &wn)
                              { surf += w * integrand(xmin + x * (xmax - xmin)); });
        // scale appropriately
        volume *= pow(xmax(0) - xmin(0), N);
        surf *= pow(xmax(0) - xmin(0), N - 1);
    };

    // Compute results for all q and output in a convergence table
    int q = qMax;
    compute(q);
    // std::cout << q << ' ' << volume << ' ' << surf << ' ' << std::abs(volume - volume_exact) / volume_exact << ' ' << std::abs(surf - surf_exact) / surf_exact << std::endl;
}
/// perform h refinement study
template <int N>
void GetFaceQuadSchemeX(uvector<real, N> xp, uvector<real, N> ymin, uvector<real, N> ymax, const uvector<int, N> &P, int q,
                        std::vector<uvector<real, N + 1>> &phase0)
{
    /// ===================== functions with type double ===================== ///
    /// calculate the level-set function value
    /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
    /// \param[out] phi - level-set function value
    auto fphi = [&](const uvector<real, 1> &xs)
    {
        using std::exp;
        using std::sqrt;
        real delta = 1e-10;
        real xscale = 1.0;
        real yscale = 1.0;
        real min_x = 0.0;
        real min_y = 0.0;
        int lsign = 1.0;
        /// dimension
        const int Ndim = N + 1;
        std::vector<uvector<real, Ndim>> xbnd;
        std::vector<uvector<real, Ndim>> norm;
        std::vector<uvector<real, Ndim - 1>> kappa;
        int nbnd = 64;
        std::cout << "nbnd " << nbnd << std::endl;
        /// radius
        real a = 2.0;
        /// center coords
        real xc = 0.0;
        real yc = 0.0;
        for (int k = 0; k < nbnd; ++k)
        {
            double theta = k * 2.0 * M_PI / nbnd;
            uvector<real, Ndim> x, nrm;
            x(0) = a * cos(theta) + xc;
            x(1) = a * sin(theta) + yc;
            nrm(0) = 2.0 * (x(0) - xc);
            nrm(1) = 2.0 * (x(1) - yc);
            real ds = nrm(0) * nrm(0) + nrm(1) * nrm(1);
            ds = sqrt(ds);
            uvector<real, Ndim> ni;
            ni(0) = nrm(0) / ds;
            ni(1) = nrm(1) / ds;
            xbnd.push_back(x);
            norm.push_back(ni);
            kappa.push_back(0.0);
        }
        uvector<real, Ndim> x;
        int rho = 10 * nbnd;
        x(0) = xs(0) * xscale + min_x;
        x(1) = xp(0) * yscale + min_y;
        /// find minimum distance
        real min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            x_diff(0) = x(0) - xbnd.at(i)(0);
            x_diff(1) = x(1) - xbnd.at(i)(1);
            real deli_x = sqrt(sqrnorm(x_diff) + delta);
            min_dist = std::min(min_dist, deli_x);
        }
        real denom = 0.0;
        real phi = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            uvector<real, Ndim> norvx, norvy;
            uvector<real, Ndim> proj;
            x_diff = x - xbnd.at(i);
            real dist1 = dot(x_diff, norm.at(i));
            real nx = norm.at(i)(0);
            real ny = norm.at(i)(1);
            // norvx = {1.0 - (nx * nx), -nx * ny};
            norvx(0) = 1.0 - (nx * nx);
            norvx(1) = -nx * ny;
            // norvy = {-nx * ny, 1.0 - (ny * ny)};
            norvy(0) = -nx * ny;
            norvy(1) = 1.0 - (ny * ny);
            proj(0) = dot(norvx, x_diff);
            proj(1) = dot(norvy, x_diff);
            real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
            real dist = dist1 + dist2;
            real delx = sqrt(sqrnorm(x_diff) + delta);
            real expc = exp(-rho * (delx - min_dist));
            denom += expc;
            phi += dist * expc;
        }
        phi = phi / denom;
        std::cout << "phi " << lsign * phi << std::endl;
        return lsign * phi;
    };
    // auto fphi = [&](const uvector<real, 1> &x)
    // {
    //     return xp(0) * xp(0) + x(0) * x(0) - 4.0;
    // };
    auto airfoil_phi = [&](const uvector<real, 1> &xs)
    {
        using std::exp;
        using std::sqrt;
        real delta = 1e-10;
        real xscale = 1.0;
        real yscale = 1.0;
        real min_x = 0.0;
        real min_y = 0.0;
        int lsign = 1.0;
        /// dimension
        const int Ndim = 2;
        std::vector<uvector<real, Ndim>> Xc;
        std::vector<uvector<real, Ndim>> xbnd;
        std::vector<uvector<real, Ndim>> norm;
        std::vector<uvector<real, Ndim - 1>> kappa;
        const int npts = 200;
        const int nbnd = 2 * npts - 3;
        // std::cout << "nbnd " << nbnd << std::endl;
        /// calculate boundary coordinates
        double tc = 0.12;
        uvector<real, npts - 1> beta;
        double beta_max = M_PI / 1.02;
        double dbeta = beta_max / (npts - 2);
        beta(0) = beta_max;
        for (int i = 1; i < npts - 1; ++i)
        {
            beta(i) = beta(i - 1) - dbeta;
        }
        constexpr bool cte = true;
        uvector<real, nbnd> xb;
        uvector<real, nbnd> yb;
        real a0 = 0.2969;
        real a1 = -0.126;
        real a2 = -0.3516;
        real a3 = 0.2843;
        real a4 = -0.1015;
        if (cte)
        {
            a4 = -0.1036;
        }
        /// upper boundary
        /// upper boundary
        for (int i = 0; i < npts; ++i)
        {
            xb(i) = (1.0 - cos(beta(i))) / 2.0;
            real term1 =
                (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
            real term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
            yb(i) = 5.0 * tc * (term1 + term2);
        }
        /// lower boundary
        for (int i = 0; i < npts - 2; ++i)
        {
            xb(i + npts - 1) = xb(npts - 3 - i);
            yb(i + npts - 1) = -yb(npts - 3 - i);
        }
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x;
            x(0) = xb(i);
            x(1) = yb(i);
            Xc.push_back(x);
        }
        /// construct the normal vector for all boundary points
        norm = constructNormal<Ndim>(Xc);
        kappa = getCurvature<Ndim, cte>(Xc);
        uvector<real, Ndim> xcent;
        xcent(0) = 19.5;
        xcent(1) = 20.0;
        for (int k = 0; k < nbnd; ++k)
        {
            uvector<real, Ndim> xtrans;
            for (int d = 0; d < Ndim; ++d)
            {
                xtrans(d) = Xc.at(k)(d) + xcent(d);
            }
            xbnd.push_back(xtrans);
        }
        // std::cout << "Airfoil coordinates: " << std::endl;
        // for (int k = 0; k < nbnd; ++k)
        // {
        //     std::cout << xbnd.at(k) << std::endl;
        // }
        /// parameters
        int ratio = 10;
        real rho = ratio * nbnd;
        uvector<real, Ndim> x;
        x(0) = xs(0) * xscale + min_x;
        x(1) = xp(0) * yscale + min_y;
        /// find minimum distance
        real min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            x_diff(0) = x(0) - xbnd.at(i)(0);
            x_diff(1) = x(1) - xbnd.at(i)(1);
            real deli_x = sqrt(sqrnorm(x_diff) + delta);
            min_dist = std::min(min_dist, deli_x);
        }
        real denom = 0.0;
        real phi = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            uvector<real, Ndim> norvx, norvy;
            uvector<real, Ndim> proj;
            x_diff = x - xbnd.at(i);
            real dist1 = dot(x_diff, norm.at(i));
            real nx = norm.at(i)(0);
            real ny = norm.at(i)(1);
            // norvx = {1.0 - (nx * nx), -nx * ny};
            norvx(0) = 1.0 - (nx * nx);
            norvx(1) = -nx * ny;
            // norvy = {-nx * ny, 1.0 - (ny * ny)};
            norvy(0) = -nx * ny;
            norvy(1) = 1.0 - (ny * ny);
            proj(0) = dot(norvx, x_diff);
            proj(1) = dot(norvy, x_diff);
            real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
            real dist = dist1 + dist2;
            real delx = sqrt(sqrnorm(x_diff) + delta);
            real expc = exp(-rho * (delx - min_dist));
            denom += expc;
            phi += dist * expc;
        }
        phi = phi / denom;
        // std::cout << "phi " << lsign * phi << std::endl;
        return lsign * phi;
    };
    // Construct phi by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phi(nullptr, P);
    algoim_spark_alloc(real, phi);
    // bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
    //                                    { return fphi(ymin + x * (ymax - ymin)); },
    //                                    phi);
    bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
                                       { return airfoil_phi(ymin + x * (ymax - ymin)); },
                                       phi);
    std::cout << "building quadrature hierarchy " << std::endl;
    // Build quadrature hierarchy
    ImplicitPolyQuadrature<N> ipquad(phi);
    // ipquad.type = OuterAggregate;
    // Compute quadrature scheme and record the nodes & weights; phase0 corresponds to
    // {phi < 0}, phase1 corresponds to {phi > 0}, and surf corresponds to {phi == 0}.
    ipquad.integrate(AlwaysGL, q, [&](const uvector<real, N> &x, real w)
                     {
                         std::cout << "q " << q << " x " << x << std::endl;
                         if (bernstein::evalBernsteinPoly(phi, x) < 0)
                             phase0.push_back(add_component(x, N, w));
                         // else
                         //     phase1.push_back(add_component(x, N, w));
                     });
}
template <int N>
void GetFaceQuadSchemeY(uvector<real, N> xp, uvector<real, N> ymin, uvector<real, N> ymax, const uvector<int, N> &P, int q,
                        std::vector<uvector<real, N + 1>> &phase0)
{
    /// ===================== functions with type double ===================== ///
    /// calculate the level-set function value
    /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
    /// \param[out] phi - level-set function value
    auto fphi = [&](const uvector<real, 1> &xs)
    {
        using std::exp;
        using std::sqrt;
        real delta = 1e-10;
        real xscale = 1.0;
        real yscale = 1.0;
        real min_x = 0.0;
        real min_y = 0.0;
        int lsign = 1.0;
        /// dimension
        const int Ndim = N + 1;
        std::vector<uvector<real, Ndim>> xbnd;
        std::vector<uvector<real, Ndim>> norm;
        std::vector<uvector<real, Ndim - 1>> kappa;
        int nbnd = 64;
        std::cout << "nbnd " << nbnd << std::endl;
        /// radius
        real a = 2.0;
        /// center coords
        real xc = 0.0;
        real yc = 0.0;
        for (int k = 0; k < nbnd; ++k)
        {
            double theta = k * 2.0 * M_PI / nbnd;
            uvector<real, Ndim> x, nrm;
            x(0) = a * cos(theta) + xc;
            x(1) = a * sin(theta) + yc;
            nrm(0) = 2.0 * (x(0) - xc);
            nrm(1) = 2.0 * (x(1) - yc);
            real ds = nrm(0) * nrm(0) + nrm(1) * nrm(1);
            ds = sqrt(ds);
            uvector<real, Ndim> ni;
            ni(0) = nrm(0) / ds;
            ni(1) = nrm(1) / ds;
            xbnd.push_back(x);
            norm.push_back(ni);
            kappa.push_back(0.0);
        }
        uvector<real, Ndim> x;
        int rho = 10 * nbnd;
        x(0) = xs(0) * xscale + min_x;
        x(1) = xp(0) * yscale + min_y;
        /// find minimum distance
        real min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            x_diff(0) = x(0) - xbnd.at(i)(0);
            x_diff(1) = x(1) - xbnd.at(i)(1);
            real deli_x = sqrt(sqrnorm(x_diff) + delta);
            min_dist = std::min(min_dist, deli_x);
        }
        real denom = 0.0;
        real phi = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            uvector<real, Ndim> norvx, norvy;
            uvector<real, Ndim> proj;
            x_diff = x - xbnd.at(i);
            real dist1 = dot(x_diff, norm.at(i));
            real nx = norm.at(i)(0);
            real ny = norm.at(i)(1);
            // norvx = {1.0 - (nx * nx), -nx * ny};
            norvx(0) = 1.0 - (nx * nx);
            norvx(1) = -nx * ny;
            // norvy = {-nx * ny, 1.0 - (ny * ny)};
            norvy(0) = -nx * ny;
            norvy(1) = 1.0 - (ny * ny);
            proj(0) = dot(norvx, x_diff);
            proj(1) = dot(norvy, x_diff);
            real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
            real dist = dist1 + dist2;
            real delx = sqrt(sqrnorm(x_diff) + delta);
            real expc = exp(-rho * (delx - min_dist));
            denom += expc;
            phi += dist * expc;
        }
        phi = phi / denom;
        std::cout << "phi " << lsign * phi << std::endl;
        return lsign * phi;
    };
    // auto fphi = [&](const uvector<real, 1> &x)
    // {
    //     return xp(0) * xp(0) + x(0) * x(0) - 4.0;
    // };
    auto airfoil_phi = [&](const uvector<real, 1> &xs)
    {
        using std::exp;
        using std::sqrt;
        real delta = 1e-10;
        real xscale = 1.0;
        real yscale = 1.0;
        real min_x = 0.0;
        real min_y = 0.0;
        int lsign = 1.0;
        /// dimension
        const int Ndim = N + 1;
        std::vector<uvector<real, Ndim>> Xc;
        std::vector<uvector<real, Ndim>> xbnd;
        std::vector<uvector<real, Ndim>> norm;
        std::vector<uvector<real, Ndim - 1>> kappa;
        const int npts = 200;
        const int nbnd = 2 * npts - 3;
        // std::cout << "nbnd " << nbnd << std::endl;
        /// calculate boundary coordinates
        double tc = 0.12;
        uvector<real, npts - 1> beta;
        double beta_max = M_PI / 1.02;
        double dbeta = beta_max / (npts - 2);
        beta(0) = beta_max;
        for (int i = 1; i < npts - 1; ++i)
        {
            beta(i) = beta(i - 1) - dbeta;
        }
        constexpr bool cte = true;
        uvector<real, nbnd> xb;
        uvector<real, nbnd> yb;
        real a0 = 0.2969;
        real a1 = -0.126;
        real a2 = -0.3516;
        real a3 = 0.2843;
        real a4 = -0.1015;
        if (cte)
        {
            a4 = -0.1036;
        }
        /// upper boundary
        /// upper boundary
        for (int i = 0; i < npts; ++i)
        {
            xb(i) = (1.0 - cos(beta(i))) / 2.0;
            real term1 =
                (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
            real term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
            yb(i) = 5.0 * tc * (term1 + term2);
        }
        /// lower boundary
        for (int i = 0; i < npts - 2; ++i)
        {
            xb(i + npts - 1) = xb(npts - 3 - i);
            yb(i + npts - 1) = -yb(npts - 3 - i);
        }
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x;
            x(0) = xb(i);
            x(1) = yb(i);
            Xc.push_back(x);
        }
        /// construct the normal vector for all boundary points
        norm = constructNormal<Ndim>(Xc);
        kappa = getCurvature<Ndim, cte>(Xc);
        uvector<real, Ndim> xcent;
        xcent(0) = 19.5;
        xcent(1) = 20.0;
        for (int k = 0; k < nbnd; ++k)
        {
            uvector<real, Ndim> xtrans;
            for (int d = 0; d < Ndim; ++d)
            {
                xtrans(d) = Xc.at(k)(d) + xcent(d);
            }
            xbnd.push_back(xtrans);
        }
        // std::cout << "Airfoil coordinates: " << std::endl;
        // for (int k = 0; k < nbnd; ++k)
        // {
        //     std::cout << xbnd.at(k) << std::endl;
        // }
        /// parameters
        int ratio = 10;
        real rho = ratio * nbnd;
        uvector<real, Ndim> x;
        x(0) = xp(0) * xscale + min_x;
        x(1) = xs(0) * yscale + min_y;
        /// find minimum distance
        real min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            x_diff(0) = x(0) - xbnd.at(i)(0);
            x_diff(1) = x(1) - xbnd.at(i)(1);
            real deli_x = sqrt(sqrnorm(x_diff) + delta);
            min_dist = std::min(min_dist, deli_x);
        }
        real denom = 0.0;
        real phi = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            uvector<real, Ndim> x_diff;
            uvector<real, Ndim> norvx, norvy;
            uvector<real, Ndim> proj;
            x_diff = x - xbnd.at(i);
            real dist1 = dot(x_diff, norm.at(i));
            real nx = norm.at(i)(0);
            real ny = norm.at(i)(1);
            // norvx = {1.0 - (nx * nx), -nx * ny};
            norvx(0) = 1.0 - (nx * nx);
            norvx(1) = -nx * ny;
            // norvy = {-nx * ny, 1.0 - (ny * ny)};
            norvy(0) = -nx * ny;
            norvy(1) = 1.0 - (ny * ny);
            proj(0) = dot(norvx, x_diff);
            proj(1) = dot(norvy, x_diff);
            real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
            real dist = dist1 + dist2;
            real delx = sqrt(sqrnorm(x_diff) + delta);
            real expc = exp(-rho * (delx - min_dist));
            denom += expc;
            phi += dist * expc;
        }
        phi = phi / denom;
        // std::cout << "phi " << lsign * phi << std::endl;
        return lsign * phi;
    };
    // Construct phi by mapping [0,1] onto bounding box [xmin,xmax]
    xarray<real, N> phi(nullptr, P);
    algoim_spark_alloc(real, phi);
    // bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
    //                                    { return fphi(ymin + x * (ymax - ymin)); },
    //                                    phi);
    bernstein::bernsteinInterpolate<N>([&](const uvector<real, N> &x)
                                       { return airfoil_phi(ymin + x * (ymax - ymin)); },
                                       phi);
    std::cout << "building quadrature hierarchy " << std::endl;
    // Build quadrature hierarchy
    ImplicitPolyQuadrature<N> ipquad(phi);
    // ipquad.type = OuterAggregate;
    // Compute quadrature scheme and record the nodes & weights; phase0 corresponds to
    // {phi < 0}, phase1 corresponds to {phi > 0}, and surf corresponds to {phi == 0}.
    ipquad.integrate(AlwaysGL, q, [&](const uvector<real, N> &x, real w)
                     {
                         std::cout << "q " << q << " x " << x << std::endl;
                         if (bernstein::evalBernsteinPoly(phi, x) < 0)
                             phase0.push_back(add_component(x, N, w));
                         // else
                         //     phase1.push_back(add_component(x, N, w));
                     });
}
#if ALGOIM_EXAMPLES_DRIVER == 0 || ALGOIM_EXAMPLES_DRIVER == 4
int main(int argc, char *argv[])
{
#if 0
    // q-convergence study for a 2D ellipse
    {
        auto ellipse = [](const uvector<real, 2> &x)
        {
            return x(0) * x(0) + x(1) * x(1) * 4 - 1;
        };
        auto integrand = [](const uvector<real, 2> &x)
        {
            return 1.0;
        };
        real volume_exact = algoim::util::pi / 2;
        real surf_exact = 4.844224110273838099214251598195914705976959198943300412541558176231060;
        std::cout << "\n\nEllipse q-convergence test\n";
        std::cout << "q      area(q)         perim(q)        area error       perim error\n";
        //
        int nel = 4;
        std::cout << "Area of a 2D ellipse, computed via the cells of a " << nel << " by " << nel << " Cartesian grid:\n";
        real dx = 2.2 / nel;
        const int N = 2;
        std::vector<uvector<real, N + 1>> surf_all, phase0_all, phase1_all;
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
                xmin(0) = -1.1 + i * dx;
                xmin(1) = -1.1 + j * dx;
                xmax(0) = -1.1 + i * dx + dx;
                xmax(1) = -1.1 + j * dx + dx;
                GetQuadScheme<2>(ellipse, xmin, xmax, 3, 10, surf, phase0, phase1);
                qAreaPeri<2>(ellipse, xmin, xmax, 3, integrand, 10, volume_exact, surf_exact, area_c, peri_c);
                // std::cout << "area: " << area_c << std::endl;
                // std::cout << "perimeter : " << peri_c << std::endl;
                area += area_c;
                peri += peri_c;
                // std::cout << "surf.size()  " << surf.size() << std::endl;
                for (int nq = 0; nq < surf.size(); ++nq)
                {
                    real xq = surf.at(nq)(0);
                    real yq = surf.at(nq)(1);
                    surf.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    surf.at(nq)(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    surf_all.push_back(surf.at(nq));
                }
                for (int nq = 0; nq < phase0.size(); ++nq)
                {
                    real xq = phase0.at(nq)(0);
                    real yq = phase0.at(nq)(1);
                    phase0.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    phase0.at(nq)(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    phase0_all.push_back(phase0.at(nq));
                }
                for (int nq = 0; nq < phase1.size(); ++nq)
                {
                    phase1_all.push_back(phase1.at(nq));
                }
            }
        std::cout << "ellipse area: " << area << std::endl;
        std::cout << "ellipse perimeter: " << peri << std::endl;
        outputQuadratureRuleAsVtpXML<N>(phase0_all, "ellipse-phase0.vtp");
        outputQuadratureRuleAsVtpXML<N>(phase1_all, "ellipse-phase1.vtp");
        outputQuadratureRuleAsVtpXML<N>(surf_all, "ellipse-surf.vtp");
    }
#endif
#if 0
     // q-convergence study for a 2D circle
    {
        auto circle = [](const uvector<real, 2> &x)
        {
            return x(0) * x(0) + x(1) * x(1) - 4.0;
        };
        auto integrand  = [](const uvector<real, 2> &x)
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
        const int N = 2;
        std::vector<uvector<real, N + 1>> surf_all, phase0_all, phase1_all;
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
                GetQuadScheme<2>(circle, xmin, xmax, 3, 10, surf, phase0, phase1);
                qAreaPeri<2>(circle, xmin, xmax, 3, integrand, 10, volume_exact, surf_exact, area_c, peri_c);
                // std::cout << "area: " << area_c << std::endl;
                // std::cout << "perimeter : " << peri_c << std::endl;
                area += area_c;
                peri += peri_c;
                // std::cout << "surf.size()  " << surf.size() << std::endl;
                for (int nq = 0; nq < surf.size(); ++nq)
                {
                    real xq = surf.at(nq)(0);
                    real yq = surf.at(nq)(1);
                    surf.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    surf.at(nq)(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    surf_all.push_back(surf.at(nq));
                }
                for (int nq = 0; nq < phase0.size(); ++nq)
                {
                    real xq = phase0.at(nq)(0);
                    real yq = phase0.at(nq)(1);
                    phase0.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    phase0.at(nq)(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    phase0_all.push_back(phase0.at(nq));
                }
                for (int nq = 0; nq < phase1.size(); ++nq)
                {
                    phase1_all.push_back(phase1.at(nq));
                }
            }
        std::cout << "circle area: " << area << std::endl;
        std::cout << "circle perimeter: " << peri << std::endl;
        outputQuadratureRuleAsVtpXML<N>(phase0_all, "circle-phase0.vtp");
        outputQuadratureRuleAsVtpXML<N>(phase1_all, "circle-phase1.vtp");
        outputQuadratureRuleAsVtpXML<N>(surf_all, "circle-surf.vtp");
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
    // q-convergence study for an airfoil
    {
        auto airfoil_phi = [](const uvector<real, 2> &xs)
        {
            using std::exp;
            using std::sqrt;
            real delta = 1e-10;
            real xscale = 1.0;
            real yscale = 1.0;
            real min_x = 0.0;
            real min_y = 0.0;
            int lsign = 1.0;
            /// dimension
            const int N = 2;
            std::vector<uvector<real, N>> Xc;
            std::vector<uvector<real, N>> xbnd;
            std::vector<uvector<real, N>> norm;
            std::vector<uvector<real, N - 1>> kappa;
            const int npts = 200;
            const int nbnd = 2 * npts - 3;
            // std::cout << "nbnd " << nbnd << std::endl;
            /// calculate boundary coordinates
            double tc = 0.12;
            uvector<real, npts - 1> beta;
            double beta_max = M_PI / 1.02;
            double dbeta = beta_max / (npts - 2);
            beta(0) = beta_max;
            for (int i = 1; i < npts - 1; ++i)
            {
                beta(i) = beta(i - 1) - dbeta;
            }
            constexpr bool cte = true;
            uvector<real, nbnd> xb;
            uvector<real, nbnd> yb;
            real a0 = 0.2969;
            real a1 = -0.126;
            real a2 = -0.3516;
            real a3 = 0.2843;
            real a4 = -0.1015;
            if (cte)
            {
                a4 = -0.1036;
            }
            /// upper boundary
            /// upper boundary
            for (int i = 0; i < npts; ++i)
            {
                xb(i) = (1.0 - cos(beta(i))) / 2.0;
                real term1 =
                    (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
                real term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
                yb(i) = 5.0 * tc * (term1 + term2);
            }
            /// lower boundary
            for (int i = 0; i < npts - 2; ++i)
            {
                xb(i + npts - 1) = xb(npts - 3 - i);
                yb(i + npts - 1) = -yb(npts - 3 - i);
            }
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x;
                x(0) = xb(i);
                x(1) = yb(i);
                Xc.push_back(x);
            }
            /// construct the normal vector for all boundary points
            norm = constructNormal<N>(Xc);
            kappa = getCurvature<N, cte>(Xc);
            uvector<real, N> xcent;
            xcent(0) = 19.5;
            xcent(1) = 20.0;
            // xcent(0) = 0.0;
            // xcent(1) = 1.0;
            for (int k = 0; k < nbnd; ++k)
            {
                uvector<real, N> xs;
                for (int d = 0; d < N; ++d)
                {
                    xs(d) = Xc.at(k)(d) + xcent(d);
                }
                xbnd.push_back(xs);
            }
            int ratio = 10;
            real rho = ratio * nbnd;
            uvector<real, N> x;
            x(0) = xs(0) * xscale + min_x;
            x(1) = xs(1) * yscale + min_y;
            /// find minimum distance
            real min_dist = 1e+100;
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                x_diff(0) = x(0) - xbnd.at(i)(0);
                x_diff(1) = x(1) - xbnd.at(i)(1);
                real deli_x = sqrt(sqrnorm(x_diff) + delta);
                min_dist = std::min(min_dist, deli_x);
            }
            real denom = 0.0;
            real phi = 0.0;
            /// evaluate the level-set
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                uvector<real, N> norvx, norvy;
                uvector<real, N> proj;
                x_diff = x - xbnd.at(i);
                real dist1 = dot(x_diff, norm.at(i));
                real nx = norm.at(i)(0);
                real ny = norm.at(i)(1);
                // norvx = {1.0 - (nx * nx), -nx * ny};
                norvx(0) = 1.0 - (nx * nx);
                norvx(1) = -nx * ny;
                // norvy = {-nx * ny, 1.0 - (ny * ny)};
                norvy(0) = -nx * ny;
                norvy(1) = 1.0 - (ny * ny);
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                real dist = dist1 + dist2;
                real delx = sqrt(sqrnorm(x_diff) + delta);
                real expc = exp(-rho * (delx - min_dist));
                denom += expc;
                phi += dist * expc;
            }
            phi = phi / denom;
            return lsign * phi;
        };
        auto integrand = [](const uvector<real, 1> &x)
        {
            return 1.0;
        };
        real volume_exact = 0.081705151533582;
        real surf_exact = 2.03955;
        std::cout << "\n\nCircle q-convergence test\n";
        std::cout << "q      area(q)         perim(q)        area error       perim error\n";
        //
        int nel = 32;
        std::cout << "Area of an airfoil, computed via the cells of a " << nel << " by " << nel << " Cartesian grid:\n";
        real dx = 1.2 / nel;
        real dy = 0.2 / nel;
        const int N = 1;
        std::vector<uvector<real, N + 2>> phaseX_all, phaseY_all;
        std::vector<uvector<real, N + 2>> surf_all, phase0_all, phase1_all;
        real area = 0.0;
        real peri = 0.0;
        for (int i = 0; i < nel; ++i)
            for (int j = 0; j < nel; ++j)
            {
                std::vector<uvector<real, N + 1>> phaseX, phaseY;
                real area_c;
                real peri_c;
                uvector<real, N + 1> xmin;
                uvector<real, N + 1> xmax;
                xmin(0) = 19.4 + i * dx;
                xmin(1) = 19.9 + j * dy;
                xmax(0) = 19.4 + i * dx + dx;
                xmax(1) = 19.9 + j * dy + dy;
                std::vector<uvector<real, N + 2>> surf, phase0, phase1;
                GetQuadScheme<2>(airfoil_phi, xmin, xmax, 11, 10, surf, phase0, phase1);
                GetFaceQuadSchemeY<1>(xmin(0), xmin(1), xmax(1), 12, 10,
                                      phaseY);
                GetFaceQuadSchemeX<1>(xmin(1), xmin(0), xmax(0), 12, 10,
                                      phaseX);
                std::vector<uvector<real, N + 2>> interface_quadY;
                for (int nq = 0; nq < surf.size(); ++nq)
                {
                    real xq = surf.at(nq)(0);
                    real yq = surf.at(nq)(1);
                    surf.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    surf.at(nq)(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    surf_all.push_back(surf.at(nq));
                }
                for (int nq = 0; nq < phaseY.size(); ++nq)
                {
                    real yq = phaseY.at(nq)(0);
                    uvector<real, N + 2> face_quad;
                    face_quad(0) = xmin(0);
                    face_quad(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    face_quad(2) = phaseY.at(nq)(1);
                    interface_quadY.push_back(face_quad);
                    phaseY_all.push_back(interface_quadY.at(nq));
                }
                std::vector<uvector<real, N + 2>> interface_quadX;
                for (int nq = 0; nq < phaseX.size(); ++nq)
                {
                    real xq = phaseX.at(nq)(0);
                    uvector<real, N + 2> face_quad;
                    face_quad(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    face_quad(1) = xmin(1);
                    face_quad(2) = phaseX.at(nq)(1);
                    interface_quadX.push_back(face_quad);
                    phaseX_all.push_back(interface_quadX.at(nq));
                }
            }
        std::cout << "airfoil area: " << area << std::endl;
        std::cout << "airfoil perimeter: " << peri << std::endl;
        outputQuadratureRuleAsVtpXML<N + 1>(phaseX_all, "airfoil-phaseX.vtp");
        outputQuadratureRuleAsVtpXML<N + 1>(phaseY_all, "airfoil-phaseY.vtp");
        outputQuadratureRuleAsVtpXML<N + 1>(surf_all, "airfoil-surf.vtp");
    }
#if 0
    std::cout << "Algoim Examples - High-order quadrature algorithms for multi-component domains implicitly-defined\n";
    std::cout << "by (one or more) multivariate Bernstein polynomials\n\n";
    std::cout << std::scientific << std::setprecision(10);
    {
        auto airfoil_phi = [](const uvector<real, 2> &xs)
        {
            using std::exp;
            using std::sqrt;
            real delta = 1e-10;
            real xscale = 1.0;
            real yscale = 1.0;
            real min_x = 0.0;
            real min_y = 0.0;
            int lsign = 1.0;
            /// dimension
            const int N = 2;
            std::vector<uvector<real, N>> Xc;
            std::vector<uvector<real, N>> xbnd;
            std::vector<uvector<real, N>> norm;
            std::vector<uvector<real, N - 1>> kappa;
            const int npts = 200;
            const int nbnd = 2 * npts - 3;
            // std::cout << "nbnd " << nbnd << std::endl;
            /// calculate boundary coordinates
            double tc = 0.12;
            uvector<real, npts - 1> beta;
            double beta_max = M_PI / 1.02;
            double dbeta = beta_max / (npts - 2);
            beta(0) = beta_max;
            for (int i = 1; i < npts - 1; ++i)
            {
                beta(i) = beta(i - 1) - dbeta;
            }
            constexpr bool cte = true;
            uvector<real, nbnd> xb;
            uvector<real, nbnd> yb;
            real a0 = 0.2969;
            real a1 = -0.126;
            real a2 = -0.3516;
            real a3 = 0.2843;
            real a4 = -0.1015;
            if (cte)
            {
                a4 = -0.1036;
            }
            /// upper boundary
            /// upper boundary
            for (int i = 0; i < npts; ++i)
            {
                xb(i) = (1.0 - cos(beta(i))) / 2.0;
                real term1 =
                    (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
                real term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
                yb(i) = 5.0 * tc * (term1 + term2);
            }
            /// lower boundary
            for (int i = 0; i < npts - 2; ++i)
            {
                xb(i + npts - 1) = xb(npts - 3 - i);
                yb(i + npts - 1) = -yb(npts - 3 - i);
            }
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x;
                x(0) = xb(i);
                x(1) = yb(i);
                Xc.push_back(x);
            }
            /// construct the normal vector for all boundary points
            norm = constructNormal<N>(Xc);
            kappa = getCurvature<N, cte>(Xc);
            uvector<real, N> xcent;
            xcent(0) = 19.5;
            xcent(1) = 20.0;
            // xcent(0) = 0.0;
            // xcent(1) = 1.0;
            for (int k = 0; k < nbnd; ++k)
            {
                uvector<real, N> xs;
                for (int d = 0; d < N; ++d)
                {
                    xs(d) = Xc.at(k)(d) + xcent(d);
                }
                xbnd.push_back(xs);
            }
            // std::cout << "Airfoil coordinates: " << std::endl;
            // for (int k = 0; k < nbnd; ++k)
            // {
            //     std::cout << xbnd.at(k) << std::endl;
            // }
            /// get the number of boundary points
            // nbnd = Xc.size();
            // std::cout << "nbnd " << nbnd << std::endl;
            // std::cout << "xbnd.size() " << xbnd.size() << std::endl;
            /// parameters
            int ratio = 10;
            real rho = ratio * nbnd;
            // std::cout << "rho " << rho << std::endl;
            // std::cout << "Normal vector: " << std::endl;
            // for (int k = 0; k < nbnd; ++k)
            // {
            //     std::cout << norm.at(k) << std::endl;
            // }
            // std::cout << "Curvature: " << std::endl;
            // for (int k = 0; k < nbnd; ++k)
            // {
            //     std::cout << kappa.at(k) << std::endl;
            // }
            uvector<real, N> x;
            x(0) = xs(0) * xscale + min_x;
            x(1) = xs(1) * yscale + min_y;
            /// find minimum distance
            real min_dist = 1e+100;
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                x_diff(0) = x(0) - xbnd.at(i)(0);
                x_diff(1) = x(1) - xbnd.at(i)(1);
                real deli_x = sqrt(sqrnorm(x_diff) + delta);
                min_dist = std::min(min_dist, deli_x);
            }
            real denom = 0.0;
            real phi = 0.0;
            /// evaluate the level-set
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                uvector<real, N> norvx, norvy;
                uvector<real, N> proj;
                x_diff = x - xbnd.at(i);
                real dist1 = dot(x_diff, norm.at(i));
                real nx = norm.at(i)(0);
                real ny = norm.at(i)(1);
                // norvx = {1.0 - (nx * nx), -nx * ny};
                norvx(0) = 1.0 - (nx * nx);
                norvx(1) = -nx * ny;
                // norvy = {-nx * ny, 1.0 - (ny * ny)};
                norvy(0) = -nx * ny;
                norvy(1) = 1.0 - (ny * ny);
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                real dist = dist1 + dist2;
                real delx = sqrt(sqrnorm(x_diff) + delta);
                real expc = exp(-rho * (delx - min_dist));
                denom += expc;
                phi += dist * expc;
            }
            phi = phi / denom;
            // std::cout << "phi " << lsign * phi << std::endl;
            return lsign * phi;
        };
        auto integrand = [](const uvector<real, 2> &x)
        {
            return 1.0;
        };
        //
        int nel = 4;
        std::cout << "Area of a 2D ellipse, computed via the cells of a " << nel << " by " << nel << " Cartesian grid:\n";
        real dx = 1.1 / nel;
        const int N = 2;
        std::vector<uvector<real, N + 1>> surf_all, phase0_all, phase1_all;
        std::vector<uvector<real, N + 1>> face0_all, face1_all;
        real area = 0.0;
        real peri = 0.0;
        real volume_exact = 0.081705151533582;
        real surf_exact = 2.03955;
        for (int i = 0; i < nel; ++i)
            for (int j = 0; j < nel; ++j)
            {
                std::vector<uvector<real, N + 1>> surf, phase0, phase1;
                std::vector<uvector<real, N + 1>> face0, face1;
                real area_c;
                real peri_c;
                uvector<real, N> xmin;
                uvector<real, N> xmax;
                xmin(0) = 19.45 + i * dx;
                xmin(1) = 19.45 + j * dx;
                xmax(0) = 19.45 + i * dx + dx;
                xmax(1) = 19.45 + j * dx + dx;
                GetQuadScheme<2>(airfoil_phi, xmin, xmax, 11, 10, surf, phase0, phase1);
                GetFaceQuadScheme<1>(airfoil_phi, xmin(0), xmax(0), 11, 10, face0, face1);
                qAreaPeri<2>(airfoil_phi, xmin, xmax, 11, integrand, 10, volume_exact, surf_exact, area_c, peri_c);
                area += area_c;
                peri += peri_c;
                // std::cout << "surf.size()  " << surf.size() << std::endl;
                for (int nq = 0; nq < surf.size(); ++nq)
                {
                    real xq = surf.at(nq)(0);
                    real yq = surf.at(nq)(1);
                    // surf.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    // surf.at(nq)(1) = xmin(1) + yq * (xmax(1) - xmin(1));
                    surf_all.push_back(surf.at(nq));
                }
                for (int nq = 0; nq < face0.size(); ++nq)
                {
                    real xq = face0.at(nq)(0);
                    real yq = face0.at(nq)(1);
                    face0.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    face0_all.push_back(face0.at(nq));
                }
                for (int nq = 0; nq < face1.size(); ++nq)
                {
                    real xq = face1.at(nq)(0);
                    real yq = face1.at(nq)(1);
                    face1.at(nq)(0) = xmin(0) + xq * (xmax(0) - xmin(0));
                    face1_all.push_back(face1.at(nq));
                }
            }
        std::cout << "Airfoil area: " << area << std::endl;
        std::cout << "Airfoil perimeter: " << peri << std::endl;
        outputQuadratureRuleAsVtpXML<N>(surf_all, "airfoilLSF.vtp");
        outputQuadratureRuleAsVtpXML<N>(face0_all, "airfoilsurf0LSF.vtp");
        outputQuadratureRuleAsVtpXML<N>(face1_all, "airfoilsurf1LSF.vtp");
    }
#endif
#if 0
    {
        auto airfoil_phiu = [](const uvector<real, 2> &xs)
        {
            using std::exp;
            using std::sqrt;
            real delta = 1e-10;
            real xscale = 1.0;
            real yscale = 1.0;
            real min_x = 0.0;
            real min_y = 0.0;
            int lsign = 1.0;
            /// dimension
            const int N = 2;
            std::vector<uvector<real, N>> Xc;
            std::vector<uvector<real, N>> xbnd;
            std::vector<uvector<real, N>> norm;
            std::vector<uvector<real, N - 1>> kappa;
            const int npts = 100;
            const int nbnd = npts;
            std::cout << "nbnd " << nbnd << std::endl;
            /// calculate boundary coordinates
            double tc = 0.12;
            uvector<real, npts - 1> beta;
            double beta_max = M_PI / 1.02;
            double dbeta = beta_max / (npts - 2);
            beta(0) = beta_max;
            for (int i = 1; i < npts - 1; ++i)
            {
                beta(i) = beta(i - 1) - dbeta;
            }
            constexpr bool cte = true;
            uvector<real, nbnd> xb;
            uvector<real, nbnd> yb;
            real a0 = 0.2969;
            real a1 = -0.126;
            real a2 = -0.3516;
            real a3 = 0.2843;
            real a4 = -0.1015;
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
                real term1 = (a0 * pow(xb(i), 0.5)) + (a1 * xb(i)) + (a2 * (pow(xb(i), 2)));
                real term2 = (a3 * (pow(xb(i), 3))) + (a4 * (pow(xb(i), 4)));
                yb(i) = 5.0 * tc * (term1 + term2);
            }
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x;
                x(0) = xb(i);
                x(1) = yb(i);
                Xc.push_back(x);
            }
            /// construct the normal vector for all boundary points
            norm = constructNormal<N>(Xc);
            kappa = getCurvature<N, cte>(Xc);
            uvector<real, N> xcent;
            xcent(0) = 19.5;
            xcent(1) = 20.0;
            // xcent(0) = 0.0;
            // xcent(1) = 0.2;
            for (int k = 0; k < nbnd; ++k)
            {
                uvector<real, N> xs;
                for (int d = 0; d < N; ++d)
                {
                    xs(d) = Xc.at(k)(d) + xcent(d);
                }
                xbnd.push_back(xs);
            }
            /// get the number of boundary points
            std::cout << "nbnd " << nbnd << std::endl;
            std::cout << "xbnd.size() " << xbnd.size() << std::endl;
            /// parameters
            int ratio = 10;
            real rho = ratio * nbnd;
            std::cout << "rho " << rho << std::endl;
            uvector<real, N> x;
            x(0) = xs(0) * xscale + min_x;
            x(1) = xs(1) * yscale + min_y;
            /// find minimum distance
            real min_dist = 1e+100;
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                x_diff(0) = x(0) - xbnd.at(i)(0);
                x_diff(1) = x(1) - xbnd.at(i)(1);
                real deli_x = sqrt(sqrnorm(x_diff) + delta);
                min_dist = std::min(min_dist, deli_x);
            }
            real denom = 0.0;
            real phi = 0.0;
            /// evaluate the level-set
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                uvector<real, N> norvx, norvy;
                uvector<real, N> proj;
                x_diff = x - xbnd.at(i);
                real dist1 = dot(x_diff, norm.at(i));
                real nx = norm.at(i)(0);
                real ny = norm.at(i)(1);
                // norvx = {1.0 - (nx * nx), -nx * ny};
                norvx(0) = 1.0 - (nx * nx);
                norvx(1) = -nx * ny;
                // norvy = {-nx * ny, 1.0 - (ny * ny)};
                norvy(0) = -nx * ny;
                norvy(1) = 1.0 - (ny * ny);
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                real dist = dist1 + dist2;
                real delx = sqrt(sqrnorm(x_diff) + delta);
                real expc = exp(-rho * (delx - min_dist));
                denom += expc;
                phi += dist * expc;
            }
            phi = phi / denom;
            std::cout << "phi " << lsign * phi << std::endl;
            return lsign * phi;
        };
        outputQuadScheme<2>(airfoil_phiu, 19.5 , 20.5, 5, 10, "airfoilUpLSF");
        // auto integrand = [](const uvector<real, 2> &x)
        // {
        //     return 1.0;
        // };
        // int qmax = 20;
        // real volume_exact = 0.0817073 / 2.0;
        // real surf_exact = 2.03955 / 2.0;
        // std::cout << "\n\nEllipse q-convergence test\n";
        // std::cout << "q      area(q)         perim(q)        area error       perim error\n";
        // qConv<2>(airfoil_phi, 0.0, 1.0, 10, integrand, qmax, volume_exact, surf_exact);
    }
#endif
/// circle
#if 0
    {
        /// ===================== functions with type double ===================== ///
        /// calculate the level-set function value
        /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
        /// \param[out] phi - level-set function value
        auto phi = [](const uvector<real, 2> &xs)
        {
            using std::exp;
            using std::sqrt;
            real delta = 1e-10;
            real xscale = 1.0;
            real yscale = 1.0;
            real min_x = 0.0;
            real min_y = 0.0;
            int lsign = 1.0;
            /// dimension
            const int N = 2;
            std::vector<uvector<real, N>> xbnd;
            std::vector<uvector<real, N>> norm;
            std::vector<uvector<real, N - 1>> kappa;
            int nbnd = 64;
            std::cout << "nbnd " << nbnd << std::endl;
            /// radius
            real a = 1.0;
            /// center coords
            real xc = 4.0;
            real yc = 4.0;
            for (int k = 0; k < nbnd; ++k)
            {
                double theta = k * 2.0 * M_PI / nbnd;
                uvector<real, N> x, nrm;
                x(0) = a * cos(theta) + xc;
                x(1) = a * sin(theta) + yc;
                nrm(0) = 2.0 * (x(0) - xc);
                nrm(1) = 2.0 * (x(1) - yc);
                real ds = nrm(0) * nrm(0) + nrm(1) * nrm(1);
                ds = sqrt(ds);
                uvector<real, N> ni;
                ni(0) = nrm(0) / ds;
                ni(1) = nrm(1) / ds;
                xbnd.push_back(x);
                norm.push_back(ni);
                kappa.push_back(0.0);
            }
            uvector<real, N> x;
            int rho = 10 * nbnd;
            x(0) = xs(0) * xscale + min_x;
            x(1) = xs(1) * yscale + min_y;
            /// find minimum distance
            real min_dist = 1e+100;
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                x_diff(0) = x(0) - xbnd.at(i)(0);
                x_diff(1) = x(1) - xbnd.at(i)(1);
                real deli_x = sqrt(sqrnorm(x_diff) + delta);
                min_dist = std::min(min_dist, deli_x);
            }
            real denom = 0.0;
            real phi = 0.0;
            /// evaluate the level-set
            for (int i = 0; i < nbnd; ++i)
            {
                uvector<real, N> x_diff;
                uvector<real, N> norvx, norvy;
                uvector<real, N> proj;
                x_diff = x - xbnd.at(i);
                real dist1 = dot(x_diff, norm.at(i));
                real nx = norm.at(i)(0);
                real ny = norm.at(i)(1);
                // norvx = {1.0 - (nx * nx), -nx * ny};
                norvx(0) = 1.0 - (nx * nx);
                norvx(1) = -nx * ny;
                // norvy = {-nx * ny, 1.0 - (ny * ny)};
                norvy(0) = -nx * ny;
                norvy(1) = 1.0 - (ny * ny);
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                real dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                real dist = dist1 + dist2;
                real delx = sqrt(sqrnorm(x_diff) + delta);
                real expc = exp(-rho * (delx - min_dist));
                denom += expc;
                phi += dist * expc;
            }
            phi = phi / denom;
            std::cout << "phi " << lsign * phi << std::endl;
            return lsign * phi;
        };
        outputQuadScheme<2>(phi, 3, 5, 3, 10, "circleLSF");
        auto integrand = [](const uvector<real, 2> &x)
        {
            return 1.0;
        };
        int qmax = 20;
        real volume_exact = M_PI;
        real surf_exact = 2.0 * M_PI;
        // std::cout << "\n\nCircle q-convergence test\n";
        // std::cout << "q      area(q)         perim(q)        area error       perim error\n";
        // qConv<2>(phi, 2.9, 5.1, 10, integrand, qmax, volume_exact, surf_exact);
    }
#endif
    return 0;
}
#endif
