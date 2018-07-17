// ======================================================================
// PureCLIP: capturing target-specific protein-RNA interaction footprints
// ======================================================================
// Copyright (C) 2017  Sabrina Krakau, Max Planck Institute for Molecular 
// Genetics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// =======================================================================
// Author: Sabrina Krakau <krakau@molgen.mpg.de>
// =======================================================================



#ifndef APPS_HMMS_DENSITY_FUNCTIONS_H_
#define APPS_HMMS_DENSITY_FUNCTIONS_H_
   
#include <iostream>
#include <fstream>
#include <math.h>       // lgamma 

#include <boost/math/tools/minima.hpp>      // BRENT's algorithm
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>       // normalized lower incomplete gamma function: gamma_p()
#include <boost/math/distributions/binomial.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

using namespace seqan;



/////////
// GAMMA2: left threshold, forced to be zero
/////////


template<typename TDOUBLE>
class GAMMA2  // ignore positions with KDE below theshold
{
public:
    GAMMA2(double theta_, double k_, double tp_): theta(theta_), k(k_), tp(tp_) {}
    GAMMA2(double tp_): tp(tp_) {}
    GAMMA2() {}

    long double getDensity(double const &x);
    void updateTheta(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const& options); 
    void updateK(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, double &kMin, double &kMax, AppOptions const& options);
    bool updateThetaAndK(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, double &kMin, double &kMax, AppOptions const& options); 

    double theta;   // scale parameter
    double mean;
    double k;       // shape parameter 
    double tp;      // truncation point
};


// Functor for Brent's algorithm: find theta
template<typename TDOUBLE>
struct FctLL_GAMMA2_theta
{
    FctLL_GAMMA2_theta(double const& k_, String<String<String<TDOUBLE> > > const& statePosteriors_,  String<String<Observations> > &setObs_, AppOptions const&options_) : k(k_), statePosteriors(statePosteriors_), setObs(setObs_), options(options_)
    { 
    }
    double operator()(double const& theta)
    {
        // Group log-likelihood function evaluations regarding binned kde vaues ! TODO 
        // normalized lower incomplete gamma function
        double nligf = boost::math::gamma_p(k, (options.useKdeThreshold/theta));
        
        TDOUBLE ll = 0.0;       
        for (unsigned s = 0; s < 2; ++s)
        {
            String<TDOUBLE> lls;
            resize(lls, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)  
                {
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1)
                    {
                        double kde = setObs[s][i].kdes[t];
                        TDOUBLE p = (k-1.0)*log(kde) - kde/theta - k*log(theta) - lgamma(k);
                        p -= log(1.0 - nligf);       
                        lls[i] +=  p * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
                ll += lls[i];
        }
        return (-ll);
    }

private:
    double k;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > &setObs;
    AppOptions options;
};


template<typename TDOUBLE>
void GAMMA2<TDOUBLE>::updateTheta(String<String<String<TDOUBLE> > > &statePosteriors, 
                         String<String<Observations> > &setObs, 
                         AppOptions const&options)
{ 
    int bits = 60;
    boost::uintmax_t maxIter = options.maxIter_brent;
    
    double thetaMin = 0.0;
    double thetaMax = 10.0;

    FctLL_GAMMA2_theta<TDOUBLE> fct_GAMMA2_theta(this->k, statePosteriors, setObs, options);
    std::pair<double, double> res = boost::math::tools::brent_find_minima(fct_GAMMA2_theta, thetaMin, thetaMax, bits, maxIter);         // use somehow initial guess to save time? or interval around prev. value?

    this->theta = res.first;
}


// Functor for Brent's algorithm: find k
template<typename TDOUBLE>
struct FctLL_GAMMA2_k
{
    FctLL_GAMMA2_k(double const& theta_, String<String<String<TDOUBLE> > > const& statePosteriors_,  String<String<Observations> > &setObs_, AppOptions const&options_) : theta(theta_), statePosteriors(statePosteriors_), setObs(setObs_), options(options_)
    { 
    }
    double operator()(double const& k)
    {
        // Group log-likelihood function evaluations regarding binned kde vaues ! TODO 
        // normalized lower incomplete gamma function
        double nligf = boost::math::gamma_p(k, (options.useKdeThreshold/theta));
 
        TDOUBLE ll = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<TDOUBLE> lls;
            resize(lls, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)  
                {
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1)
                    {
                        double kde = setObs[s][i].kdes[t];
                        TDOUBLE p = (k-1.0)*log(kde) - kde/theta - k*log(theta) - lgamma(k);
                        p -= log(1.0 - nligf);
                        lls[i] +=  p * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
                ll += lls[i];
        } 
        return (-ll);
    }

private:
    double theta;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > &setObs;
    AppOptions options;
};


template<typename TDOUBLE>
void GAMMA2<TDOUBLE>::updateK(String<String<String<TDOUBLE> > > &statePosteriors, 
                     String<String<Observations> > &setObs,
                     double &kMin, double &kMax,
                     AppOptions const&options)
{ 
    int bits = 60;
    boost::uintmax_t maxIter = options.maxIter_brent;
    
    FctLL_GAMMA2_k<TDOUBLE> fct_GAMMA2_k(this->theta, statePosteriors, setObs, options);
    std::pair<double, double> res = boost::math::tools::brent_find_minima(fct_GAMMA2_k, kMin, kMax, bits, maxIter);         // use somehow initial guess to save time? or interval around prev. value?

    this->k = res.first;
}



/// use GSL simplex to update k and theta together
template<typename TDOUBLE>
struct Fct_GSL_X_GAMMA2
{
    Fct_GSL_X_GAMMA2(double const & tp_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_,  String<String<Observations> > &setObs_, 
                                  AppOptions const&options_) : tp(tp_), 
                                                               statePosteriors(statePosteriors_), 
                                                               setObs(setObs_), 
                                                               options(options_)
    { 
    }
    // f
    double operator()(const gsl_vector * x)
    {       
        const double theta = gsl_vector_get (x, 0);
        const double k = gsl_vector_get (x, 1);
        double pred = k*theta;

        if (k <= 0) std::cerr << "ERROR calling GSL multimin: k = " << k << " (make sure initial step size is small enough, value not getting <= 0!)" << std::endl;
        if (theta <= 0) std::cerr << "ERROR calling GSL multimin: theta = " << theta << " (make sure initial step size is small enough, value not getting <= 0!)" << std::endl;

        double nligf = boost::math::gamma_p(k, (tp*theta));

        TDOUBLE f = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<TDOUBLE> f_S;
            resize(f_S, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)
                {    
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1) 
                    {
                        double kde = setObs[s][i].kdes[t];
            
                        TDOUBLE p = (k-1.0)*log(kde) - k * (kde/pred + log(pred)) - k*log(1.0/k) - lgamma(k) - log(1.0 - nligf);

                        f_S[i] +=  p * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
                f += f_S[i];
        }
        return  (-f);  
    }
   
private:
    double tp;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > & setObs;
    AppOptions options;
};


template<typename TDOUBLE>
struct Fct_GSL_X_GAMMA2_fixK
{
    Fct_GSL_X_GAMMA2_fixK(double const & tp_, double const& k_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_,  String<String<Observations> > &setObs_, 
                                  AppOptions const&options_) : tp(tp_), k(k_),
                                                               statePosteriors(statePosteriors_),  
                                                               setObs(setObs_),  
                                                               options(options_)
    { 
    }
    // f
    double operator()(const gsl_vector * x)
    {       
        const double theta = gsl_vector_get (x, 0);

        double pred = k*theta;

        if (k <= 0) std::cerr << "ERROR calling GSL multimin: k = " << k << " (make sure initial step size is small enough, value not getting <= 0!)" << std::endl;
        if (theta <= 0) std::cerr << "ERROR calling GSL multimin: theta = " << theta << " (make sure initial step size is small enough, value not getting <= 0!)" << std::endl;

        double nligf = boost::math::gamma_p(k, (tp*theta));

        TDOUBLE f = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<double> f_S;
            resize(f_S, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)
                {    
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1) 
                    {
                        double kde = setObs[s][i].kdes[t];
           
                        TDOUBLE p = (k-1.0)*log(kde) - k * (kde/pred + log(pred)) - k*log(1.0/k) - lgamma(k) - log(1.0 - nligf);

                        f_S[i] +=  p * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
                f += f_S[i];
        }
        return  (-f);  
    }
   
private:
    double tp;
    double k;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > &setObs;
    AppOptions options;
};

// Wrapper functions for functors
template<typename TDOUBLE>
double fct_GSL_X_GAMMA2_W (const gsl_vector * x, void * p) {

    Fct_GSL_X_GAMMA2<TDOUBLE> * function = reinterpret_cast< Fct_GSL_X_GAMMA2<TDOUBLE> *> (p);
    return (*function)( x );        
} 

template<typename TDOUBLE>
double fct_GSL_X_GAMMA2_fixK_W (const gsl_vector * x, void * p) {

    Fct_GSL_X_GAMMA2_fixK<TDOUBLE> * function = reinterpret_cast< Fct_GSL_X_GAMMA2_fixK<TDOUBLE> *> (p);
    return (*function)( x );        
} 


template<typename TDOUBLE>
struct Params3
{
    double tp;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > setObs;
    AppOptions options;
};


template<typename TDOUBLE>
struct Params4
{
    double tp;
    double k;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > setObs;
    AppOptions options;
};


template<typename TDOUBLE>
bool callGSL_simplex2_fixK(int &status, 
                  double &tp, double &theta, double &k,
                  String<String<String<TDOUBLE> > > &statePosteriors, 
                  String<String<Observations> > &setObs, 
                  AppOptions const& options)
{
    int iter = 0;
    int max_iter = 100;
    const size_t n = 1; 
    double size;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    
    struct Params4<TDOUBLE> params = {tp, k, statePosteriors, setObs, options};
    gsl_multimin_function f;

    // instantiation of functor with all fixed params
    Fct_GSL_X_GAMMA2_fixK<TDOUBLE> fct(tp, k, statePosteriors, setObs, options);

    /* Set initial step sizes to */
    gsl_vector *ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 0.001);  


    f.n = n;
    f.f = &fct_GSL_X_GAMMA2_fixK_W<TDOUBLE>;        // pointer to wrapper member function
    f.params =  &fct;       // pointer to functor (instead of to params)

    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, theta);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, n);
    gsl_multimin_fminimizer_set (s, &f, x, ss);  

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);  

        if (status)
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-6);

        if (options.verbosity >= 2)
        {
            if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

            printf ("%5d %10.7f f() = %7.7f size = %.7f\n", 
                  iter,
                  gsl_vector_get (s->x, 0), 
                  s->fval, size);
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    theta = gsl_vector_get (s->x, 0);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free(ss);
    return true;
}


template<typename TDOUBLE>
bool callGSL_simplex2(double &tp, double &theta, double &k,
                  String<String<String<TDOUBLE> > > &statePosteriors, 
                  String<String<Observations> > &setObs, 
                  double &kMin, double &kMax,
                  AppOptions const& options)
{
    if (options.verbosity >= 2) 
        std::cout << "Call GSL multimin solver nmsimplex2 ..." << std::endl;

    int status;
    int iter = 0;
    int max_iter = 200;
    const size_t n = 2; 
    double size;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    
    struct Params3<TDOUBLE> params = {tp, statePosteriors, setObs, options};
    gsl_multimin_function f;

    // instantiation of functor with all fixed params
    Fct_GSL_X_GAMMA2<TDOUBLE> fct(tp, statePosteriors, setObs, options);

    /* Set initial step sizes to 0.0001 */
    gsl_vector *ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 0.001);  // TODO different for k and theta (make sure not getting below 0!)
    // TODO adjust to given value 

    f.n = n;
    f.f = &fct_GSL_X_GAMMA2_W<TDOUBLE>;        // pointer to wrapper member function
    f.params =  &fct;       // pointer to functor (instead of to params)
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, theta);
    gsl_vector_set (x, 1, k);
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, n);
    gsl_multimin_fminimizer_set (s, &f, x, ss);   
    
    if (options.verbosity >= 2)
    {
        printf ("%5d %10.7f %10.7f", 
                  0,
                  gsl_vector_get (s->x, 0), 
                  gsl_vector_get (s->x, 1));
    }

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);  

        if (status)
        break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-6);

        if (options.verbosity >= 2)
        {
            if (status == GSL_SUCCESS)
                printf ("Minimum found at:\n");
        }
        // if k < kMin: fix k and optimize only for theta
        if (gsl_vector_get (s->x, 1) < kMin)
        {
            std::cout << "Note: limited shape parameter k to: " << kMin << ". Make sure g1.k <= g2.k. Decrease in shape could be caused by outliers: high peaks, potentially background binding regions. Check if transcripts/chromosomes used for learning are representative. Incorporating input signal would help. (usually results still show relatively high precision compared to other methods)" <<  std::endl;

            theta = gsl_vector_get (s->x, 0);
            callGSL_simplex2_fixK(status, tp, theta, kMin, statePosteriors, setObs, options);  

            gsl_vector_set (s->x, 0, theta);
            gsl_vector_set (s->x, 1, kMin);
            break;
        }
        else if (gsl_vector_get (s->x, 1) > kMax)
        {
            std::cout << "Note: limited shape parameter k to: " << kMax << std::endl;
            theta = gsl_vector_get (s->x, 0);
            callGSL_simplex2_fixK(status, tp, theta, kMax, statePosteriors, setObs, options);  

            gsl_vector_set (s->x, 0, theta);
            gsl_vector_set (s->x, 1, kMax);
            break;
        }

        if (options.verbosity >= 2)
        {
            printf ("%5d %10.7f %10.7f f() = %7.7f size = %.7f\n", 
                  iter,
                  gsl_vector_get (s->x, 0), 
                  gsl_vector_get (s->x, 1), 
                  s->fval, size);
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if (options.verbosity >= 2)
    {
        printf ("status = %s\n", gsl_strerror (status));
        std::cout << "GSL simplex2 .... theta = " << gsl_vector_get (s->x, 0)  << " k = " << gsl_vector_get (s->x, 1) << std::endl;
    }
    theta = gsl_vector_get (s->x, 0);
    k = gsl_vector_get (s->x, 1); 

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free(ss);
    return true;
}



template<typename TDOUBLE>
bool GAMMA2<TDOUBLE>::updateThetaAndK(String<String<String<TDOUBLE> > > &statePosteriors, 
                    String<String<Observations> > &setObs, 
                    double &kMin, double &kMax,
                    AppOptions const&options)
{
    // use multidimensional simplex2
    return callGSL_simplex2(this->tp, this->theta, this->k, statePosteriors, setObs, kMin, kMax, options);    
}



//
///////////////////////////////////////////////


template<typename TDOUBLE>
long double GAMMA2<TDOUBLE>::getDensity(double const &x)   
{
    if (x < this->tp) return 0.0; 

    TDOUBLE f1 = pow(x, this->k - 1.0) * exp(-x/this->theta);
    TDOUBLE f2 = pow(this->theta, this->k) * tgamma(this->k);
    
    // normalized lower incomplete gamma function
    double nligf = boost::math::gamma_p(this->k, this->tp/this->theta);

    return  ( (f1/f2) / (1.0 - nligf));
}



//////////////////////////
// utils


template<typename TDOUBLE>
bool loadGammaParams(GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2, AppOptions &options)
{
    std::ifstream ifs(toCString(options.inParFileName));
    if (!ifs.good())
        std::cerr << "ERROR: Could not open file containing input parameters!\n";
    std::string lineBuffer;

    while (std::getline(ifs, lineBuffer))
    {
        //std::cout << lineBuffer << std::endl;
        std::string value1;
        std::string value2;
        std::stringstream ss(lineBuffer);
        if (!ss.str().empty())
        {
            if (!std::getline(ss, value1, '\t'))
            {
                std::cerr << "ERROR: could not read first value.\n";
                return false;
            }
            if (value1.c_str() == "gamma1.theta")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for gamma1.theta\n";
                    return false;
                }
                gamma1.theta = std::strtod(value2.c_str(), NULL);
            }
            else if (value1.c_str() == "gamma1.k")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for gamma1.k\n";
                    return false;
                }
                gamma1.k = std::strtod(value2.c_str(), NULL);
            }
            else if (value1.c_str() == "gamma1.tp")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for gamma1.tp\n";
                    return false;
                }
                gamma1.tp = std::strtod(value2.c_str(), NULL);
            }
            else if (value1.c_str() == "gamma2.theta")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for gamma2.theta\n";
                    return false;
                }
                gamma2.theta = std::strtod(value2.c_str(), NULL);
            }
            else if (value1.c_str() == "gamma2.k")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for gamma2.k\n";
                    return false;
                }
                gamma2.k = std::strtod(value2.c_str(), NULL);
            }
            else if (value1.c_str() == "gamma2.tp")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for gamma2.tp\n";
                    return false;
                }
                gamma2.tp = std::strtod(value2.c_str(), NULL);
            }
        }
    }
    return true;
}



template<typename TDOUBLE>
void myPrint(GAMMA2<TDOUBLE> &gamma)
{
    std::cout << "*** GAMMA2 ***" << std::endl;
    std::cout << "    theta:"<< gamma.theta << std::endl;
    std::cout << "    k:" << gamma.k << std::endl;
    std::cout << "    tp:" << gamma.tp << std::endl;
    std::cout << std::endl;
}

template<typename TDOUBLE>
bool checkConvergence(GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2, AppOptions &options)
{
    if (std::fabs(gamma1.theta - gamma2.theta) > options.gamma_theta_conv) return false;
    if (std::fabs(gamma1.k - gamma2.k) > options.gamma_k_conv) return false;

    return true;
}

template<typename TOut, typename TDOUBLE>
void printParams(TOut &out, GAMMA2<TDOUBLE> &gamma, int i)
{
    out << "gamma" << i << ".theta" << '\t' << gamma.theta << std::endl;
    out << "gamma" << i << ".k" << '\t' << gamma.k << std::endl;
    out << "gamma" << i << ".tp" << '\t' << gamma.tp << std::endl;
    out << std::endl;
}


template<typename TDOUBLE>
void checkOrderG1G2(GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2, AppOptions &options)
{
    if ((gamma1.k*gamma1.theta) > (gamma2.k*gamma2.theta))
    {
        std::swap(gamma1.theta, gamma2.theta);
        std::cout << "NOTE: swapped gamma1.theta and gamma2.theta ! " << std::endl;
    }
}


#endif
