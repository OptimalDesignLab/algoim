#ifndef ALGOIM_LEVELSET_HPP
#define ALGOIM_LEVELSET_HPP
#include <array>
#include <bitset>
#include <vector>
#include <algorithm>
#include "algoim_blitzinc.hpp"
#include "algoim_utility.hpp"
#include "algoim_quad.hpp"
#include <chrono>
using namespace std::chrono;
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
                                vector<TinyVector<double, N - 1>> _kappa, double _rho, double _lsign = -1.0, double _delta = 1e-10)
        {
            xbnd = _xbnd;
            norm = _norm;
            kappa = _kappa;
            rho = _rho;
            delta = _delta;
            lsign = _lsign;
            xscale = 1.0;
            yscale = 1.0;
            min_x = 0.0;
            min_y = 0.0;
            sign_phi = lsign;
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
                TinyVector<double, N> norvx, norvy;
                TinyVector<double, N> proj;
                x_diff = x - xbnd.at(i);
                double dist1 = dot(x_diff, norm.at(i));
                double nx = norm.at(i)(0);
                double ny = norm.at(i)(1);
                norvx = {1.0 - (nx * nx), -nx * ny};
                norvy = {-nx * ny, 1.0 - (ny * ny)};
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                double dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                double dist = dist1 + dist2;
                double delx = sqrt(magsqr(x_diff) + delta);
                double expc = exp(-rho * (delx - min_dist));
                denom += expc;
                phi += dist * expc;
            }
            phi = phi / denom;
            return lsign * phi;
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
                TinyVector<double, N> norvx, norvy;
                TinyVector<double, N> proj;
                double perp1 = dot(x_diff, norm.at(i));
                double nx = norm.at(i)(0);
                double ny = norm.at(i)(1);
                norvx = {1.0 - (nx * nx), -nx * ny};
                norvy = {-nx * ny, 1.0 - (ny * ny)};
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                double perp2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                double perp = perp1 + perp2;
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
                TinyVector<double, N> norvx, norvy;
                TinyVector<double, N> proj;
                double perp1 = dot(x_diff, norm.at(i));
                double nx = norm.at(i)(0);
                double ny = norm.at(i)(1);
                norvx = {1.0 - (nx * nx), -nx * ny};
                norvy = {-nx * ny, 1.0 - (ny * ny)};
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                double perp2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                double perp = perp1 + perp2;
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
                // add curvature correction
                phi_bar += perp_bar * kappa.at(i)(0) * proj;
                // dist = sqrt(dot(x - xc, x - xc) + levset.delta)
                phi_bar += (dist_bar / dist) * x_diff;
            }
            // cout << "phi_bar " << phi_bar << endl;
            return lsign * phi_bar;
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
                TinyVector<double, N> norvx, norvy;
                TinyVector<double, N> proj;
                x_diff = x - xbnd.at(i);
                double perp1 = dot(x_diff, norm.at(i));
                double nx = norm.at(i)(0);
                double ny = norm.at(i)(1);
                norvx = {1.0 - (nx * nx), -nx * ny};
                norvy = {-nx * ny, 1.0 - (ny * ny)};
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                double perp2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                double perp = perp1 + perp2;
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
                TinyVector<double, N> norvx, norvy;
                TinyVector<double, N> proj;
                double perp1 = dot(x_diff, norm.at(i));
                double nx = norm.at(i)(0);
                double ny = norm.at(i)(1);
                norvx = {1.0 - (nx * nx), -nx * ny};
                norvy = {-nx * ny, 1.0 - (ny * ny)};
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                double perp2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);
                double perp = perp1 + perp2;
                double delx = sqrt(magsqr(x_diff) + delta);
                double expc = exp(-rho * (delx - min_dist));
                TinyVector<double, N> perp2_xv;
                perp2_xv = kappa.at(i)(0) * proj;
                for (int d = 0; d < N; ++d)
                {
                    double expc_x = -rho * expc * x_diff(d) / delx;
                    double perp1_bar = norm.at(i)(d);
                    double perp2_bar = perp2_xv(d);
                    double perp_bar = perp1_bar + perp2_bar;
                    numer_x(d) += (perp * expc_x) + (expc * perp_bar);
                    denom_x(d) += expc_x;
                    term1 = delx * (expc_x * x_diff(d) + expc);
                    term2 = expc * (x_diff(d) * x_diff(d)) / delx;
                    double expc_xx = -rho * (term1 - term2) / (delx * delx);
                    /// this is zero for linear case
                    double perp1_xx = 0.0;
                    /// linear for quadratic distance correction
                    double nxv = norm.at(i)(d);
                    double perp2_xx = kappa.at(i)(0) * (1.0 - (nxv * nxv));
                    double perp_xx = perp1_xx + perp2_xx;
                    term1 = (perp * expc_xx) + (2.0 * perp_bar * expc_x);
                    term2 = expc * perp_xx;
                    numer_xx(d) += term1 + term2;
                    denom_xx(d) += expc_xx;
                }
                double expc_x = -rho * expc * x_diff(0) / delx;
                double expc_y = -rho * expc * x_diff(1) / delx;
                double perp1_x = norm.at(i)(0);
                double perp1_y = norm.at(i)(1);
                double perp2_x = perp2_xv(0);
                double perp2_y = perp2_xv(1);
                double perp_x = perp1_x + perp2_x;
                double perp_y = perp1_y + perp2_y;
                term1 = delx * expc_y * x_diff(0);
                term2 = expc * x_diff(0) * x_diff(1) / delx;
                double expc_xy = -rho * (term1 - term2) / (delx * delx);
                term1 = delx * expc_x * x_diff(1);
                term2 = expc * x_diff(0) * x_diff(1) / delx;
                double expc_yx = -rho * (term1 - term2) / (delx * delx);
                double perp1_xy = 0.0;
                double perp1_yx = 0.0;
                double perp2_xy = -kappa.at(i)(0) * norm.at(i)(0) * norm.at(i)(1);
                double perp2_yx = -kappa.at(i)(0) * norm.at(i)(0) * norm.at(i)(1);
                double perp_xy = perp1_xy + perp2_xy;
                double perp_yx = perp1_yx + perp2_yx;
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
        #if 1
        // template <typename T>
        Interval<N> operator()(const blitz::TinyVector<Interval<N>, N> &x) const
        {
            using std::sqrt;
            using T = Interval<N>;
            int nbnd = xbnd.size();
            TinyVector<double, N> xc;
            xc(0) = x(0).alpha;
            xc(1) = x(1).alpha;
            TinyVector<double, N> dx = x(0).delta();
            double L = sqrt(magsqr(dx));
            /// find minimum distance
            double min_dist = 1e+100;
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<double, N> x_diff;
                x_diff = xc - xbnd.at(i);
                double deli_x = sqrt(magsqr(x_diff) + delta);
                min_dist = min(min_dist, deli_x);
            }
            T denom = T(0.0);
            T phi = T(0.0);
            /// evaluate the level-set
            for (int i = 0; i < nbnd; ++i)
            {
                TinyVector<T, N> x_diff;
                TinyVector<double, N> norvx, norvy;
                TinyVector<T, N> proj;
                x_diff = x - xbnd.at(i);
                T dist1 = dot(x_diff, norm.at(i));
                double nx = norm.at(i)(0);
                double ny = norm.at(i)(1);

                /// curvature correction
                norvx = {1.0 - (nx * nx), -nx * ny};
                norvy = {-nx * ny, 1.0 - (ny * ny)};
                proj(0) = dot(norvx, x_diff);
                proj(1) = dot(norvy, x_diff);
                T dist2 = 0.5 * kappa.at(i)(0) * dot(x_diff, proj);

                T dist = dist1 + dist2;
                TinyVector<double, N> x_d;
                x_d = xc - xbnd.at(i);
                double del_xc = sqrt(magsqr(x_d) + delta);
                TinyVector<double, N> beta;
                for (int dim = 0; dim < N; ++dim)
                {
                    beta(dim) = x_d(dim)/del_xc;
                }
                double dist_max = max(sqrt(delta), del_xc - L);
                double eps = min(2.0*L, L*L/dist_max);
                T delx = Interval<N>(del_xc, beta, eps);
                T expc = exp(-rho * (delx - min_dist));
                denom += expc;
                phi += dist * expc;
            }
            //cout << "denom " << denom.alpha << endl;
            //   cout << "problem ? " << endl;
            phi = phi / denom;
            //cout << "phi " << endl;
            // cout << phi.alpha << endl;
            // cout << "phi bounds " << endl;
            cout <<  phi.maxDeviation() << endl;
            return lsign * phi;
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
        /// count 
        // static double cphi = 0.0;
        double sign_phi;

    private:
        /// Vector of boundary coordinates
        vector<TinyVector<double, N>> xbnd;
        /// Vector of boundary normal vectors
        vector<TinyVector<double, N>> norm;
        /// curvature at every point
        vector<TinyVector<double, N - 1>> kappa;
        /// penalty parameter
        double rho;
        /// sign of ls
        double lsign;
        /// parameter that smooths distance near zero
        double delta;
    }; // namespace Algoim
} // namespace Algoim
#endif