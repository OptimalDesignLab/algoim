#ifndef ALGOIM_LEVELSET_HPP
#define ALGOIM_LEVELSET_HPP
#include <array>
#include <bitset>
#include <vector>
#include <algorithm>
#include "algoim_blitzinc.hpp"
#include "algoim_utility.hpp"
#include "algoim_quad.hpp"
using namespace blitz;
using namespace std;
using namespace Algoim;
namespace Algoim
{
template <int N>
class LevelSet
{
public:
    /// Constructor for levelset
    /// \param[in] _xbnd - vector of the boundary coordinates
    /// \param[in] _norm - vector of boundary normal vectors
    /// \param[in] _rho  - parameter
    /// \param[in] _delta - parameter to prevent small/negative sqrt values
    // LevelSet(vector<TinyVector<double, N>> _xbnd, vector<TinyVector<double, N>> _norm,
    //          double _rho, double _delta = 1e-10)
    //     : xbnd(_xbnd), norm(_norm), rho(_rho), delta(_delta) {}

    void initializeLevelSet(vector<TinyVector<double, N>> _xbnd, vector<TinyVector<double, N>> _norm,
                            double _rho, double _delta = 1e-10)
    {
        xbnd = _xbnd;
        norm = _norm;
        rho = _rho;
        delta = _delta;
    }
    /// ===================== functions with type double ===================== ///
    /// calculate the level-set function value
    /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
    /// \param[out] phi - level-set function value
    double operator()(const blitz::TinyVector<double, N> &xs) const
    {
        using std::exp;
        using std::sqrt;
        TinyVector<double, N> x;
        x[0] = xs[0] * xscale + min_x;
        x[1] = xs[1] * yscale + min_y;
        int nbnd = xbnd.size();
        /// find minimum distance
        double min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double deli_x = sqrt(magsqr(x_diff) + delta);
            min_dist = min(min_dist, deli_x);
        }
        double denom = 0.0;
        double phi = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double dist = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            denom += expc;
            phi += dist * expc;
        }
        phi = phi / denom;
        return -phi;
    }

    /// calculate the gradient of level-set function
    /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
    /// \param[out] phi_bar - level-set function gradient value
    TinyVector<double, N> grad(const blitz::TinyVector<double, N> &xs) const
    {
        using std::exp;
        using std::sqrt;
        TinyVector<double, N> x;
        x[0] = xs[0] * xscale + min_x;
        x[1] = xs[1] * yscale + min_y;
        int nbnd = xbnd.size();
        /// find minimum distance
        double min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double deli_x = sqrt(magsqr(x_diff) + delta);
            min_dist = min(min_dist, deli_x);
        }
        double numer = 0.0;
        double denom = 0.0;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double perp = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            denom += expc;
            numer += perp * expc;
        }
        double ls = numer / denom;
        // start reverse sweep
        // return ls
        double ls_bar = 1.0;
        // ls = numer / denom
        double numer_bar = ls_bar / denom;
        double denom_bar = -(ls_bar * ls) / denom;
        TinyVector<double, N> phi_bar;
        phi_bar = 0.0;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double dist = sqrt(magsqr(x_diff) + delta);
            double perp = dot(x_diff, norm.at(i));
            double expfac = exp(-rho * (dist - min_dist));
            // denom += expfac
            double expfac_bar = denom_bar;
            expfac_bar += numer_bar * perp;
            // numer += perp*expfac
            double perp_bar = numer_bar * expfac;
            // expfac = exp(-rho*dist)
            double dist_bar = -expfac_bar * expfac * rho;
            // perp = dot(levset.normal[:,i], x - xc)
            phi_bar += perp_bar * norm.at(i);
            // dist = sqrt(dot(x - xc, x - xc) + levset.delta)
            phi_bar += (dist_bar / dist) * x_diff;
        }
        //cout << "phi_bar " << phi_bar << endl;
        return -phi_bar;
    }

    /// calculate the second derivative of level-set function
    /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
    /// \param[out] phi_xx - level-set function's second derivative
    TinyVector<double, 2 * N> hessian(const blitz::TinyVector<double, N> &x) const
    {
        using std::exp;
        using std::sqrt;
        // TinyVector<double, N> x;
        // x[0] = xs[0] * xscale + xmin;
        // x[1] = xs[1] * yscale + ymin;
        int nbnd = xbnd.size();
        TinyVector<double, N> phi_xx;
        TinyVector<double, 2 * N> phix_bar;
        /// find minimum distance
        double min_dist = 1e+100;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double deli_x = sqrt(magsqr(x_diff) + delta);
            min_dist = min(min_dist, deli_x);
        }
        double numer = 0.0;
        double denom = 0.0;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double perp = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            denom += expc;
            numer += perp * expc;
        }
        TinyVector<double, N> numer_x = 0.0;
        TinyVector<double, N> numer_xx = 0.0;
        TinyVector<double, N> denom_x = 0.0;
        TinyVector<double, N> denom_xx = 0.0;
        double numer_xy = 0.0;
        double denom_xy = 0.0;
        double numer_yx = 0.0;
        double denom_yx = 0.0;
        double term1;
        double term2;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = x - xbnd.at(i);
            double perp = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            for (int d = 0; d < N; ++d)
            {
                double expc_x = -rho * expc * x_diff(d) / delx;
                double perp_bar = norm.at(i)(d);
                numer_x(d) += (perp * expc_x) + (expc * perp_bar);
                denom_x(d) += expc_x;
                term1 = delx * (expc_x * x_diff(d) + expc);
                term2 = expc * (x_diff(d) * x_diff(d)) / delx;
                double expc_xx = -rho * (term1 - term2) / (delx * delx);

                /// this is zero for linear case
                double perp_xx = 0.0;
                term1 = (perp * expc_xx) + (2.0 * perp_bar * expc_x);
                term2 = expc * perp_xx;
                numer_xx(d) += term1 + term2;
                denom_xx(d) += expc_xx;
            }

            double expc_x = -rho * expc * x_diff(0) / delx;
            double expc_y = -rho * expc * x_diff(1) / delx;
            double perp_x = norm.at(i)(0);
            double perp_y = norm.at(i)(1);
            term1 = delx * expc_y * x_diff(0);
            term2 = expc * x_diff(0) * x_diff(1) / delx;
            double expc_xy = -rho * (term1 - term2) / (delx * delx);
            term1 = delx * expc_x * x_diff(1);
            term2 = expc * x_diff(0) * x_diff(1) / delx;
            double expc_yx = -rho * (term1 - term2) / (delx * delx);
            double perp_xy = 0.0;
            double perp_yx = 0.0;
            term1 = (perp * expc_xy) + (expc_x * perp_y);
            term2 = (expc_y * perp_x) + (perp_xy * expc);
            numer_xy += term1 + term2;
            denom_xy += expc_xy;
            term1 = (perp * expc_yx) + (expc_y * perp_x);
            term2 = (expc_x * perp_y) + (perp_yx * expc);
            numer_yx += term1 + term2;
            denom_yx += expc_yx;
        }
        for (int d = 0; d < N; ++d)
        {
            term1 = (denom * numer_xx(d)) - (numer * denom_xx(d));
            term2 = (denom * numer_x(d)) - (numer * denom_x(d));
            double fac = denom * denom;
            phi_xx(d) = term1 / fac - (2.0 * term2 * denom_x(d) / (fac * denom));
        }
        term1 = (denom * numer_xy) + (numer_x(0) * denom_x(1));
        term2 = (numer * denom_xy) + (denom_x(0) * numer_x(1));
        double term3 = ((denom * numer_x(0)) - (numer * denom_x(0))) * 2.0 * denom_x(1);
        double fac = denom * denom;
        double phi_xy = (term1 - term2) / fac;
        phi_xy -= term3 / (fac * denom);
        term1 = (denom * numer_yx) + (numer_x(1) * denom_x(0));
        term2 = (numer * denom_yx) + (denom_x(1) * numer_x(0));
        term3 = 2 * ((denom * numer_x(1)) - (numer * denom_x(1))) * denom_x(0);
        double phi_yx = (term1 - term2) / fac;
        phi_yx -= term3 / (fac * denom);
        phix_bar(0) = phi_xx(0);
        phix_bar(1) = phi_xy;
        phix_bar(2) = phi_yx;
        phix_bar(3) = phi_xx(1);
        return phix_bar;
    }

    /// ===================== Interval type functions start here =====================///
    /// calculate the level-set function value (here the data type is interval)
    /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
    /// \param[out] phi - level-set function value
    Interval<N> operator()(const blitz::TinyVector<Interval<N>, N> &xs) const
    {
        using std::exp;
        using std::sqrt;
        TinyVector<Interval<N>, N> x;
        x[0] = xs[0] * xscale + min_x;
        x[1] = xs[1] * yscale + min_y;
        int nbnd = xbnd.size();
        /// find minimum distance
        double min_dist = 1e+100;
        /// get the centroid
        TinyVector<double, N> xc;
        xc[0] = x(0).alpha;
        xc[1] = x(1).alpha;
        /// find minimum distance
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = xc - xbnd.at(i);
            double deli_x = sqrt(magsqr(x_diff) + delta);
            min_dist = min(min_dist, deli_x);
        }
        double denom = 0.0;
        double phi_xc = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = xc - xbnd.at(i);
            double dist = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            denom += expc;
            phi_xc += dist * expc;
        }
        phi_xc = phi_xc / denom;
        /// get the element domain
        TinyVector<double, N> xmin;
        TinyVector<double, N> xmax;
        xmax = xc + x(0).delta();
        xmin = xc - x(0).delta();
        const int nv = pow(2, N);
        /// get element vertices
        TinyVector<TinyVector<double, N>, nv> vert;
        // \Note- need to modify this
        vert(0) = xmin;
        vert(1) = {xmax(0), xmin(1)};
        vert(2) = xmax;
        vert(3) = {xmin(0), xmax(1)};
        TinyVector<double, N> dx = x(0).delta();
        double L = sqrt(magsqr(dx));
        double fac = exp(2.0 * rho * L);
        /// evaluate level-set bounds
        double phi_bnd = 0.0;
        for (int i = 0; i < nbnd; ++i)
        {
            double dmin = 1e100;
            double dmax = -1e100;
            for (int k = 0; k < nv; ++k)
            {
                TinyVector<double, N> ds;
                ds = vert(k) - xbnd.at(i);
                double perp = dot(ds, norm.at(i));
                dmin = min(dmin, perp);
                dmax = max(dmax, perp);
            }
            TinyVector<double, N> x_diff;
            x_diff = xc - xbnd.at(i);
            double dist = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            double psi_xc = expc / denom;

            double min_psi = min(psi_xc * fac, 1.0);
            double max_phi = max(dmax - phi_xc, -dmin + phi_xc);
            phi_bnd += max_phi * min_psi;
        }
        // cout << phi_bnd << endl;
        LevelSet<N> ls;
        ls.initializeLevelSet(xbnd, norm, rho, delta);
        /// need to scale back to reference element space
        TinyVector<double, N> xcs;
        xcs(0) = (xc(0) - min_x) / xscale;
        xcs(1) = (xc(1) - min_y) / yscale;
        TinyVector<double, N> beta = -ls.grad(xcs);
        double eps = phi_bnd;
        for (int dim = 0; dim < N; ++dim)
        {
            eps -= std::abs(beta(dim)) * x(dim).delta(dim);
        }
        Interval<N> phi = Interval<N>(phi_xc, beta, eps);
        return -phi;
    }

    /// calculate the gradient of level-set function
    /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
    /// \param[out] phi_bar - level-set function gradient value

    TinyVector<Interval<N>, N> grad(const blitz::TinyVector<Interval<N>, N> &xs) const
    {
        using std::exp;
        using std::sqrt;
        // scale the x values for physical space
        TinyVector<Interval<N>, N> x;
        x[0] = xs[0] * xscale + min_x;
        x[1] = xs[1] * yscale + min_y;
        TinyVector<double, N> xc;
        xc(0) = x(0).alpha;
        xc(1) = x(1).alpha;
        int nbnd = xbnd.size();
        /// find minimum distance
        double min_dist = 1e+100;
        /// find `j` index at which psi(xc) is largest
        int j_max = 0;
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = xc - xbnd.at(i);
            double deli_x = sqrt(magsqr(x_diff) + delta);
            // min_dist = min(min_dist, deli_x);
            if (deli_x < min_dist)
            {
                min_dist = deli_x;
                j_max = i;
            }
        }
        double denom = 0.0;
        /// evaluate the level-set
        for (int i = 0; i < nbnd; ++i)
        {
            TinyVector<double, N> x_diff;
            x_diff = xc - xbnd.at(i);
            double dist = dot(x_diff, norm.at(i));
            double delx = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (delx - min_dist));
            denom += expc;
        }
        /// calculate levelset derivative bounds
        TinyVector<double, N> xmin;
        TinyVector<double, N> xmax;
        xmax = xc + x(0).delta();
        xmin = xc - x(0).delta();
        const int nv = pow(2, N);

        TinyVector<TinyVector<double, N>, nv> vert;
        // \Note- need to modify this
        vert(0) = xmin;
        vert(1) = {xmax(0), xmin(1)};
        vert(2) = xmax;
        vert(3) = {xmin(0), xmax(1)};

        double djmin = 1e100;
        double djmax = -1e100;
        for (int k = 0; k < nv; ++k)
        {
            TinyVector<double, N> ds;
            ds = vert(k) - xbnd.at(j_max);
            double perp = dot(ds, norm.at(j_max));
            djmin = min(djmin, perp);
            djmax = max(djmax, perp);
        }
        TinyVector<double, N> dx = x(0).delta();
        double L = sqrt(magsqr(dx));
        LevelSet<N> ls;
        ls.initializeLevelSet(xbnd, norm, rho, delta);

        /// scale back to reference element
        TinyVector<double, N> xcs;
        xcs(0) = (xc(0) - min_x) / xscale;
        xcs(1) = (xc(1) - min_y) / yscale;
        TinyVector<double, N> beta = -ls.grad(xcs);

        double delx_phi = beta(0);
        double dely_phi = beta(1);

        double phix_bnd = 0.0;
        double phiy_bnd = 0.0;
        for (int j = 0; j < nbnd; ++j)
        {
            double delx_perp = norm.at(j)(0);
            double dely_perp = norm.at(j)(1);
            TinyVector<double, N> x_diff;
            x_diff = xc - xbnd.at(j);
            double deli_x = sqrt(magsqr(x_diff) + delta);
            double expc = exp(-rho * (deli_x - min_dist));
            double psi_xc = expc / denom;
            double bar_psi = min(1.0, psi_xc * exp(2 * rho * L));
            phix_bnd += abs(delx_perp - delx_phi) * bar_psi;
            phiy_bnd += abs(dely_perp - dely_phi) * bar_psi;
            if (j != j_max)
            {
                double dmin = 1e100;
                double dmax = -1e100;
                double psi_bar = psi_xc * exp(-2.0 * rho * L);
                double dij_max = -1e100;
                for (int k = 0; k < nv; ++k)
                {
                    TinyVector<double, N> ds;
                    ds = vert(k) - xbnd.at(j);
                    double perp = dot(ds, norm.at(j));
                    dmin = min(dmin, perp);
                    dmax = max(dmax, perp);
                    TinyVector<double, N> ds_jmax;
                    ds_jmax = vert(k) - xbnd.at(j_max);
                    double perp_jmax = dot(ds_jmax, norm.at(j_max));
                    double d_diff = abs(perp - perp_jmax);
                    dij_max = max(d_diff, dij_max);
                }
                double min_psi = min(bar_psi * (1.0 - psi_bar), 0.25);
                double temp;
                if ((bar_psi >= 0.5) && (psi_bar <= 0.5))
                {
                    temp = 2.0 * rho * 0.25;
                }
                else
                {
                    temp = 2.0 * rho * max(psi_bar * (1.0 - psi_bar), bar_psi * (1.0 - bar_psi));
                }

/// old bounds
#if (0)
                phix_bnd += max(dmax - djmin, -dmin + djmax) * temp;
                phiy_bnd += max(dmax - djmin, -dmin + djmax) * temp;
#endif
                /// new bounds
                phix_bnd += dij_max * temp;
                phiy_bnd += dij_max * temp;
            }
        }
        // cout << phix_bnd << " , " << phiy_bnd << endl;
        double eps_x = phix_bnd;
        double eps_y = phiy_bnd;
        TinyVector<double, 2 * N> hes;
        hes = hessian(xc);
        TinyVector<double, N> beta_x, beta_y;
        beta_x(0) = hes(0);
        beta_x(1) = hes(1);
        beta_y(0) = hes(2);
        beta_y(1) = hes(3);
        for (int dim = 0; dim < N; ++dim)
        {
            eps_x -= std::abs(beta_x(dim)) * x(dim).delta(dim);
            eps_y -= std::abs(beta_y(dim)) * x(dim).delta(dim);
        }
        Interval<N> phi_x = Interval<N>(delx_phi, beta_x, eps_x);
        Interval<N> phi_y = Interval<N>(dely_phi, beta_y, eps_y);
        // Interval<N> phi_x = Interval<N>(delx_phi, beta_x);
        // Interval<N> phi_y = Interval<N>(dely_phi, beta_y);
        return blitz::TinyVector<Interval<N>, N>(-phi_x, -phi_y);
    }
#if 0
        blitz::TinyVector<Interval<N>, N> grad(const blitz::TinyVector<Interval<N>, N> &x) const
        {
            using std::exp;
            using std::sqrt;
            double term1, term2;
            /// get the centroid
            TinyVector<double, N> xc;
            xc(0) = x(0).alpha;
            xc(1) = x(1).alpha;
            TinyVector<double, N> delta_x, delta_y;
            delta_x = x(0).delta();
            TinyVector<double, N> phix_xc;
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
            // cout << "xmax " << xmax << endl;
            TinyVector<double, N> phi_g, beta_x, beta_y;
            TinyVector<double, N *(N + 1) / 2> hes;
            hes = hessian(xmax);
            beta_x(0) = hes(0);
            beta_x(1) = hes(1);
            beta_y(0) = hes(1);
            beta_y(1) = hes(2);
            phi_g = grad(xmax);
            term1 = phi_g(0) * dx;
            term2 = hes(1) * dy;
            double phix_bnd = abs(term1 + term2);
            term1 = phi_g(1) * dy;
            term2 = hes(1) * dx;
            double phiy_bnd = abs(term1 + term2);
            double eps_x = phix_bnd;
            double eps_y = phiy_bnd;
            // cout << phix_bnd << " , " << phiy_bnd << endl;
            for (int dim = 0; dim < N; ++dim)
            {
                eps_x -= std::abs(beta_x(dim)) * x(dim).delta(dim);
                eps_y -= std::abs(beta_y(dim)) * x(dim).delta(dim);
            }
            // Interval<N> phi_x = Interval<N>(phix_xc(0), beta_x, eps_x);
            // Interval<N> phi_y = Interval<N>(phix_xc(1), beta_y, eps_y);
            Interval<N> phi_x = Interval<N>(phix_xc(0), beta_x);
            Interval<N> phi_y = Interval<N>(phix_xc(1), beta_y);
            return blitz::TinyVector<Interval<N>, N>(phi_x, phi_y);
        }
#endif
#if 0
        TinyVector<Interval<N>, N> grad(const blitz::TinyVector<Interval<N>, N> &x) const
        {
            int nbnd = xbnd.size();
            /// find minimum distance
            double min_dist = 0.0;
            Interval<N> expc = Interval<N>(0);
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<Interval<N>, N> x_diff;
                x_diff = x - xbnd.at(i);
                Interval<N> deli_x = sqrt(magsqr(x_diff) + delta);
                expc += exp(-rho * (deli_x - min_dist));
                //min_dist = min(min_dist, deli_x);
            }
            Interval<N> numer = Interval<N>(0);
            Interval<N> ls = Interval<N>(0);
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<Interval<N>, N> x_diff;
                x_diff = x - xbnd.at(i);
                Interval<N> perp = dot(x_diff, norm.at(i));
                Interval<N> delx = sqrt(magsqr(x_diff) + delta);
                Interval<N> logpsi = -rho * (delx - min_dist);
                logpsi -= log(expc);
                Interval<N> psi = exp(logpsi);
                ls += perp * psi;
            }
            // start reverse sweep
            // return ls
            Interval<N> ls_bar = Interval<N>(1.0);
            // ls = numer / denom
            Interval<N> lognumer_bar = log(ls_bar) - log(expc);
            Interval<N> numer_bar = exp(lognumer_bar);
            Interval<N> logdenom_bar = -log(expc);
            Interval<N> denom_bar = -ls * ls_bar * exp(logdenom_bar);
            TinyVector<Interval<N>, N> phi_bar;
            Interval<N> phi_x = Interval<N>(0);
            Interval<N> phi_y = Interval<N>(0);
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<Interval<N>, N> x_diff;
                x_diff = x - xbnd.at(i);
                Interval<N> dist = sqrt(magsqr(x_diff) + delta);
                Interval<N> perp = dot(x_diff, norm.at(i));
                Interval<N> expfac = exp(-rho * (dist - min_dist));
                // denom += expfac
                Interval<N> expfac_bar = denom_bar;
                expfac_bar += numer_bar * perp;
                // numer += perp*expfac
                Interval<N> perp_bar = numer_bar * expfac;
                // expfac = exp(-rho*dist)
                Interval<N> dist_bar = -expfac_bar * expfac * rho;
                // perp = dot(levset.normal[:,i], x - xc)
                phi_x += perp_bar * norm.at(i)(0);
                phi_y += perp_bar * norm.at(i)(1);
                // dist = sqrt(dot(x - xc, x - xc) + levset.delta)
                Interval<N> logdist = -log(dist);
                phi_x += (dist_bar * exp(logdist)) * x_diff(0);
                phi_y += (dist_bar * exp(logdist)) * x_diff(1);
                // phi_bar += (dist_bar / dist) * x_diff;
            }
            return blitz::TinyVector<Interval<N>, N>(phi_x, phi_y);
        }
#endif
    mutable double xscale;
    mutable double yscale;
    mutable double min_x;
    mutable double min_y;

private:
    /// Vector of boundary coordinates
    vector<TinyVector<double, N>> xbnd;
    /// Vector of boundary normal vectors
    vector<TinyVector<double, N>> norm;
    /// curvature at every point
    vector<TinyVector<double, N>> kappa;
    /// penalty parameter
    double rho;
    /// parameter that smooths distance near zero
    double delta;

}; // namespace Algoim
} // namespace Algoim
#endif