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
        LevelSet(vector<TinyVector<double, N>> _xbnd, vector<TinyVector<double, N>> _norm,
                 double _rho, double _delta = 1e-10)
            : xbnd(_xbnd), norm(_norm), rho(_rho), delta(_delta) {}

        /// calculate the level-set function value
        /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
        /// \param[out] phi - level-set function value
        double operator()(const blitz::TinyVector<double, N> &x) const
        {
            using std::exp;
            using std::sqrt;
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
            return phi;
        }

        /// calculate the gradient of level-set function
        /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
        /// \param[out] phi_bar - level-set function gradient value
        TinyVector<double, N> grad(const blitz::TinyVector<double, N> &x) const
        {
            using std::exp;
            using std::sqrt;
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
            return phi_bar;
        }

        /// calculate the second derivative of level-set function
        /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
        /// \param[out] phi_xx - level-set function's second derivative
        TinyVector<double, N + 1> hessian(const blitz::TinyVector<double, N> &x) const
        {
            using std::exp;
            using std::sqrt;
            int nbnd = xbnd.size();
            TinyVector<double, N> phi_xx;
            TinyVector<double, N + 1> phix_bar;
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
                    double expc_x = rho * expc * x_diff(d) / delx;
                    double perp_bar = norm.at(i)(d);
                    numer_x(d) += (perp * expc_x) + (expc * perp_bar);
                    denom_x(d) += expc_x;
                    term1 = delx * (expc_x * x_diff(d)) + expc;
                    term2 = expc * (x_diff(d) * x_diff(d)) / delx;
                    double expc_xx = rho * (term1 + term2) / (delx * delx);

                    /// this is zero for linear case
                    double perp_xx = 0.0;
                    term1 = (perp * expc_xx) + (2.0 * perp_bar * expc_x);
                    term2 = expc * perp_xx;
                    numer_xx(d) += term1 + term2;
                    denom_xx(d) += expc_xx;
                }

                double expc_x = rho * expc * x_diff(0) / delx;
                double expc_y = rho * expc * x_diff(1) / delx;
                double perp_x = norm.at(i)(0);
                double perp_y = norm.at(i)(1);
                term1 = delx * expc_y * x_diff(0);
                term2 = expc * x_diff(0) * x_diff(1) / delx;
                double expc_xy = rho * (term1 + term2) / (delx * delx);
                double perp_xy = 0.0;
                term1 = (perp * expc_xy) + (expc_x * perp_y);
                term2 = (expc_y * perp_x) + (perp_xy * expc);
                numer_xy += term1 + term2;
                denom_xy += expc_xy;
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
            phix_bar(0) = phi_xx(0);
            phix_bar(1) = phi_xy;
            phix_bar(2) = phi_xx(1);
            return phix_bar;
        }
        /// calculate the level-set function value (here the data type is interval)
        /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
        /// \param[out] phi - level-set function value
        Interval<N> operator()(const blitz::TinyVector<Interval<N>, N> &x) const
        {
            using std::exp;
            using std::sqrt;
            int nbnd = xbnd.size();
            /// find minimum distance
            double min_dist = 1e+100;
            /// get the centroid
            TinyVector<double, N> xc;
            xc(0) = x(0).alpha;
            xc(1) = x(1).alpha;
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
            LevelSet<N> ls(xbnd, norm, rho, delta);
            TinyVector<double, N> beta = ls.grad(xc);
            double eps = phi_bnd;
            for (int dim = 0; dim < N; ++dim)
            {
                eps -= std::abs(beta(dim)) * x(0).delta(dim);
            }
            Interval<N> phi = Interval<N>(phi_xc, beta, eps);
            return phi;
        }
        
        /// calculate the gradient of level-set function
        /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
        /// \param[out] phi_bar - level-set function gradient value
        TinyVector<Interval<N>, N> grad(const blitz::TinyVector<Interval<N>, N> &x) const
        {
            using std::exp;
            using std::sqrt;
            TinyVector<double, N> xc;
            xc(0) = x(0).alpha;
            xc(1) = x(1).alpha;
            int nbnd = xbnd.size();
            /// find minimum distance
            double min_dist = 1e+100;
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
            LevelSet<N> ls(xbnd, norm, rho, delta);
            TinyVector<double, N> beta = ls.grad(xc);
            double delx_phi = beta(0);
            double dely_phi = beta(1);
            /// find `j` index at which psi(xc) is largest
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
                    double psi_bar = psi_xc * exp(-2 * rho * L);
                    for (int k = 0; k < nv; ++k)
                    {
                        TinyVector<double, N> ds;
                        ds = vert(k) - xbnd.at(j);
                        double perp = dot(ds, norm.at(j));
                        dmin = min(dmin, perp);
                        dmax = max(dmax, perp);
                    }
                    double min_psi = min(bar_psi * (1 - psi_bar), 0.25);
                    double temp;
                    if ((bar_psi > 0.5) && (psi_bar < 0.5))
                    {
                        temp = 2.0 * rho * 0.25;
                    }
                    else
                    {
                        temp = 2.0 * rho * max(psi_bar * (1.0 - psi_bar), bar_psi * (1.0 - bar_psi));
                    }
                    phix_bnd += max(dmax - djmin, -dmin + djmax) * temp;
                    phiy_bnd += max(dmax - djmin, -dmin + djmax) * 2 * rho * min_psi;
                }
            }
            double eps_x = phix_bnd;
            double eps_y = phiy_bnd;
            TinyVector<double, N + 1> hes;
            hes = hessian(xc);
            TinyVector<double, N> beta_x, beta_y;
            beta_x(0) = hes(0);
            beta_x(1) = hes(1);
            beta_y(0) = hes(1);
            beta_y(1) = hes(2);
            for (int dim = 0; dim < N; ++dim)
            {
                eps_x -= std::abs(beta_x(dim)) * x(0).delta(dim);
                eps_y -= std::abs(beta_y(dim)) * x(0).delta(dim);
            }
            Interval<N> phi_x = Interval<N>(delx_phi, beta_x, eps_x);
            Interval<N> phi_y = Interval<N>(dely_phi, beta_y, eps_y);
            return blitz::TinyVector<Interval<N>, N>(phi_x, phi_y);
        }
#if 0
        /// This checks the level-set function along a given x slice
        /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
        /// \param[out] phi - level-set function value
        Interval<N> operator()(const blitz::TinyVector<Interval<N>, N> &x) const
        {
            using std::exp;
            using std::sqrt;
            int nbnd = xbnd.size();
            const int nel = 80;
            TinyVector<double, nel> phi_bnds;
            TinyVector<double, nel> phix_bnds;
            TinyVector<double, nel> phiy_bnds;
            TinyVector<double, N> dxp;
            double ymin = -2.0;
            double ymax = 2.0;
            double xslice = 3.0;
            double ds = (ymax - ymin) / nel;
            dxp = {0.5 * ds, 0.5 * ds};
            double L = mag(dxp);
            for (int k = 0; k < nel; ++k)
            {
                TinyVector<double, N> xc;
                TinyVector<double, N> yelem;
                xc(0) = xslice;
                yelem(0) = (ymax - ymin) * (k) / nel + ymin;
                yelem(1) = (ymax - ymin) * (k + 1) / nel + ymin;
                xc(1) = 0.5 * (yelem(0) + yelem(1));
                double ymid = xc(1);
                const int nv = pow(2, N);
                /// create the element vertices
                std::vector<TinyVector<double, N>> vert;
                TinyVector<double, N> xmid = {xslice, ymid};
                for (int idx = 0; idx < nv; ++idx)
                {
                    vert.push_back(xmid);
                }
                for (int idx = 0; idx < nv; ++idx)
                {
                    for (int d = 0; d < N; ++d)
                    {
                        if ((d == 0) && ((idx == 1) || (idx == nv - 1)))
                        {

                            vert.at(idx)(d) += dxp[d];
                        }
                        else if (((d == 1) && ((idx == 2) || (idx == (nv - 1)))))
                        {
                            vert.at(idx)(d) += dxp[d];
                        }
                        else
                        {

                            vert.at(idx)(d) -= dxp[d];
                        }
                    }
                }
                int j_max = 0;

                /// find minimum distance
                double min_dist = 1e+100;
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

                //cout << "min_dist " << min_dist << endl;
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
                double fac = exp(2.0 * rho * L);
                double phi_bnd = 0.0;
                for (int i = 0; i < nbnd; ++i)
                {
                    double dmin = 1e100;
                    double dmax = -1e100;
                    for (int kv = 0; kv < nv; ++kv)
                    {
                        TinyVector<double, N> ds;
                        ds = vert.at(kv) - xbnd.at(i);
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
                    // double min_psi = min(psi_xc , 1.0/fac) * fac;
                    double min_psi = min(psi_xc * fac, 1.0);
                    double max_phi = max(dmax - phi_xc, -dmin + phi_xc);
                    phi_bnd += max_phi * min_psi;
                }

                /// verify the derivative bounds
                double djmin = 1e100;
                double djmax = -1e100;
                for (int kv = 0; kv < nv; ++kv)
                {
                    TinyVector<double, N> ds;
                    ds = vert.at(kv) - xbnd.at(j_max);
                    double perp = dot(ds, norm.at(j_max));
                    djmin = min(djmin, perp);
                    djmax = max(djmax, perp);
                }
                TinyVector<double, N> beta = grad(xc);
                double delx_phi = beta(0);
                double dely_phi = beta(1);
                /// find `j` index at which psi(xc) is largest
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
                    double bar_psi = min(1.0, psi_xc * exp(2.0 * rho * L));
                    phix_bnd += abs(delx_perp - delx_phi) * bar_psi;
                    phiy_bnd += abs(dely_perp - dely_phi) * bar_psi;
                    if (j != j_max)
                    {
                        double dmin = 1e100;
                        double dmax = -1e100;
                        double psi_bar = psi_xc * exp(-2.0 * rho * L);
                        for (int k = 0; k < nv; ++k)
                        {
                            TinyVector<double, N> ds;
                            ds = vert.at(k) - xbnd.at(j);
                            double perp = dot(ds, norm.at(j));
                            dmin = min(dmin, perp);
                            dmax = max(dmax, perp);
                        }
                        double min_psi = min(bar_psi * (1 - psi_bar), 0.25);
                        double temp;
                        if ((bar_psi > 0.5) && (psi_bar < 0.5))
                        {
                            temp = 2.0 * rho * 0.25;
                        }
                        else
                        {
                            temp = 2.0 * rho * max(psi_bar * (1.0 - psi_bar), bar_psi * (1.0 - bar_psi));
                        }
                        phix_bnd += max(dmax - djmin, -dmin + djmax) * temp;
                        phiy_bnd += max(dmax - djmin, -dmin + djmax) * 2 * rho * min_psi;
                    }
                }
                phi_bnds(k) = phi_bnd;
                phix_bnds(k) = phix_bnd;
                phiy_bnds(k) = phiy_bnd;
            }
            cout << " =============================  " << endl;
            cout << "level set bounds " << endl;
            cout << " =============================  " << endl;
            cout << phi_bnds << endl;
            cout << " =============================  " << endl;
            cout << "level set gradient bounds (x dir) " << endl;
            cout << " =============================  " << endl;
            cout << phix_bnds << endl;
            double phi_xc = 0.0;
            Interval<N> phi = Interval<N>(phi_xc, 0.0, 0.0);
            return phi;
        }
#endif

#if 0
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
#endif

/// Just checking if the Interval type functions work correctly using double
#if 0
        /// calculate the level-set function value (here the data type is interval)
        /// \param[in] x - tinyvector of point where phi(x) needs to be calculated
        /// \param[out] phi - level-set function value
        double operator()(const blitz::TinyVector<double, N> &x) const
        {
            using std::exp;
            using std::log;
            using std::sqrt;
            int nbnd = xbnd.size();
            /// find minimum distance
            double min_dist = double(0);
            double expc = 0.0;
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x_diff;
                x_diff = x - xbnd.at(i);
                double deli_x = sqrt(magsqr(x_diff) + delta);
                expc += exp(-rho * (deli_x - min_dist));
                //min_dist = min(min_dist, deli_x);
            }
            double phi = double(0);
            /// evaluate the level-set
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x_diff;
                x_diff = x - xbnd.at(i);
                double dist = dot(x_diff, norm.at(i));
                double delx = sqrt(magsqr(x_diff) + delta);
                // Interval<N> expc = exp(-rho * (delx - min_dist));
                // denom += expc;
                double logpsi = -rho * (delx - min_dist);
                logpsi -= log(expc);
                double psi = exp(logpsi);
                phi += dist * psi;
            }
            return phi;
        }

        // Just checking 
        /// calculate the gradient of level-set function
        /// \param[in] x - tinyvector of point where gradphi(x) needs to be calculated
        /// \param[out] phi_bar - level-set function gradient value
        TinyVector<double, N> grad(const blitz::TinyVector<double, N> &x) const
        {
            using std::exp;
            using std::log;
            using std::sqrt;
            int nbnd = xbnd.size();
            /// find minimum distance
            double min_dist = 0.0;
            double expc = 0.0;
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x_diff;
                x_diff = x - xbnd.at(i);
                double deli_x = sqrt(magsqr(x_diff) + delta);
                expc += exp(-rho * (deli_x - min_dist));
                //min_dist = min(min_dist, deli_x);
            }
            double numer = double(0);
            double ls = double(0);
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x_diff;
                x_diff = x - xbnd.at(i);
                double perp = dot(x_diff, norm.at(i));
                double delx = sqrt(magsqr(x_diff) + delta);
                double logpsi = -rho * (delx - min_dist);
                logpsi -= log(expc);
                double psi = exp(logpsi);
                ls += perp * psi;
            }
            // start reverse sweep
            // return ls
            double ls_bar = double(1.0);
            // ls = numer / denom
            double lognumer_bar = log(ls_bar) - log(expc);
            double numer_bar = exp(lognumer_bar);
            double logdenom_bar = -log(expc);
            double denom_bar = -ls * ls_bar * exp(logdenom_bar);
            TinyVector<double, N> phi_bar;
            double phi_x = double(0);
            double phi_y = double(0);
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
                phi_x += perp_bar * norm.at(i)(0);
                phi_y += perp_bar * norm.at(i)(1);
                // dist = sqrt(dot(x - xc, x - xc) + levset.delta)
                double logdist = -log(dist);
                phi_x += (dist_bar * exp(logdist)) * x_diff(0);
                phi_y += (dist_bar * exp(logdist)) * x_diff(1);
                // phi_bar += (dist_bar / dist) * x_diff;
            }
            return blitz::TinyVector<double, N>(phi_x, phi_y);
        }
#endif

    private:
        /// Vector of boundary coordinates
        vector<TinyVector<double, N>> xbnd;
        /// Vector of boundary normal vectors
        vector<TinyVector<double, N>> norm;
        /// parameters
        double rho;
        double delta;
    }; // namespace Algoim
} // namespace Algoim
#endif