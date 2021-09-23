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
