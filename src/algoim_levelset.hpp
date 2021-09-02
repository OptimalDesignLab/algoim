#ifndef ALGOIM_LEVELSET_HPP
#define ALGOIM_LEVELSET_HPP
#include <array>
#include <bitset>
#include <vector>
#include <algorithm>
#include "algoim_blitzinc.hpp"
#include "algoim_utility.hpp"
using namespace blitz;
using namespace std;
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
            return phi_bar;
        }
    private:
        /// Vector of boundary coordinates
        vector<TinyVector<double, N>> xbnd;
        /// Vector of boundary normal vectors
        vector<TinyVector<double, N>> norm;
        /// parameters
        double rho;
        double delta;
    };
} // namespace Algoim
#endif