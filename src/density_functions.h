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

    GAMMA2(double tp_): tp(tp_) {}
    GAMMA2() {}

    long double getDensity(double const &x);
    bool updateThetaAndK(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, double &kMin, double &kMax, AppOptions const& options); 
    bool updateThetaAndK(String<String<double> > &startSet, String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, double &kMin, double &kMax, AppOptions const& options); 

    double b0;   // scale parameter
    double k;       // shape parameter 
    double tp;      // truncation point
};



///////

template<typename TDOUBLE>
long double my_GSL_X_GAMMA2_forK(const gsl_vector * x, long double const & k, 
        long double const & tp,
        String<String<String<TDOUBLE> > > const& statePosteriors,
        String<String<Observations> > & setObs,  
        AppOptions const&options)
{       
    const long double b0 = gsl_vector_get (x, 1);   //theta
    long double pred = exp(b0); //k*theta;

    long double nligf = boost::math::gamma_p(k, (tp*k/pred));    // changed ...
    if (nligf == 1.0) 
    {
        SEQAN_OMP_PRAGMA(critical)
        if (options.verbosity >= 2) std::cout << "NOTE: nligf: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << nligf << std::setprecision(6) << " gamma mean: " << pred << std::endl;
    }

    long double f = 0.0;
    for (unsigned s = 0; s < 2; ++s)
    {
        String<long double> f_S;
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
                        long double kde = setObs[s][i].kdes[t];

                        long double p = (k-1.0)*log(kde) - k * (kde/pred + log(pred)) - k*log(1.0/k) - lgamma(k) - log(1.0 - nligf);

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



/// use GSL simplex to update k and theta together
template<typename TDOUBLE>
struct Fct_GSL_X_GAMMA2
{
    Fct_GSL_X_GAMMA2(double const & tp_,
                                  double const & minK_,
                                  double const & maxK_,
                                  double & penalty_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_,  String<String<Observations> > &setObs_, 
                                  AppOptions const&options_) : tp(tp_),
                                                               minK(minK_),
                                                               maxK(maxK_),
                                                               penalty(penalty_),
                                                               statePosteriors(statePosteriors_), 
                                                               setObs(setObs_), 
                                                               options(options_)
    { 
    }
    // f
    long double operator()(const gsl_vector * x)
    {       
        const long double k = gsl_vector_get (x, 0);

        long double f = 0.0;

        if (k >= minK && k <= maxK)                                                                 // if valid k
        {
            f = my_GSL_X_GAMMA2_forK(x, k, tp, statePosteriors, setObs, options);
        }
        else if (k < minK)
        {
            //std::cout << "k < kmin " << k << std::endl;
            long double f_c = my_GSL_X_GAMMA2_forK(x, minK, tp, statePosteriors, setObs, options);                // f value at constraint
            long double f_cn = my_GSL_X_GAMMA2_forK(x, (minK+0.001), tp, statePosteriors, setObs, options);       // f value inside the constraints with distance of 0.001
            long double d = minK - k;

            // descending towards constraint:
            // -> mirror function values at constraint line - penalty
            // only if mirror point < maxK!
            if (f_cn - f_c > 0.0 && (minK + d <= maxK))
            {
                //std::cout << "k < kmin " << k << " descending towards constraint" << std::endl;
                f = my_GSL_X_GAMMA2_forK(x, (minK+d), tp, statePosteriors, setObs, options);    // NOTE: f is already negative
                f += pow(d*(-f)*penalty, 2.0);                                                      // penalty depending on distance to constraint -> prevent simplex from moving outside of constraints   
            }
            // ascending towards constraint:
            // -> use function values at constraint line - penalty
            else // if (f_cn - f_c >= 0)
            {
                //std::cout << "k < kmin " << k << " ascending towards constraint" << std::endl;
                f = my_GSL_X_GAMMA2_forK(x, minK, tp, statePosteriors, setObs, options);
                f += pow(d*(-f)*penalty, 2.0);
            }
        }
        else                                                                                                    //if (k > maxK)
        {
            long double f_c = my_GSL_X_GAMMA2_forK(x, maxK, tp, statePosteriors, setObs, options);                // f value at constraint
            long double f_cn = my_GSL_X_GAMMA2_forK(x, (maxK-0.001), tp, statePosteriors, setObs, options);       // f value inside the constraints with distance of 0.001
            long double d = k - maxK;

            // descending towards constraint:
            // -> mirror function values at constraint line - penalty
            // only if mirror point > minK!
            if (f_cn - f_c > 0.0 && (maxK - d >= minK))
            {
                f = my_GSL_X_GAMMA2_forK(x, (maxK-d), tp, statePosteriors, setObs, options);
                f += pow(d*(-f)*penalty, 2.0);
            }
            // ascending towards constraint:
            // -> use function values at constraint line - penalty
            else // if (f_cn - f_c >= 0)
            {
                f = my_GSL_X_GAMMA2_forK(x, maxK, tp, statePosteriors, setObs, options); 
                f += pow(d*(-f)*penalty, 2.0);
            } 
        }
        return  f;  
    }
   
private:
    long double tp;
    long double minK;
    long double maxK;
    long double penalty;
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
    long double operator()(const gsl_vector * x)
    {       
        const long double b0 = gsl_vector_get (x, 0);   // theta

        long double pred = exp(b0);     // k*theta;

        long double nligf = boost::math::gamma_p((long double)k, (long double)tp*k/pred);       // ...changed
        if (nligf == 1.0) 
        {
            SEQAN_OMP_PRAGMA(critical)
            if (options.verbosity >= 2) std::cout << "NOTE: nligf: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << nligf << std::setprecision(6) << " gamma mean: " << pred << std::endl;
        }

        TDOUBLE f = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<long double> f_S;
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
                        long double kde = setObs[s][i].kdes[t];
           
                        long double p = (k-1.0)*log(kde) - k * (kde/pred + log(pred)) - k*log(1.0/k) - lgamma(k) - log(1.0 - nligf);

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
    long double tp;
    long double k;
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
    double minK;
    double maxK;
    double penalty;
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
                  double &fval,
                  double &tp, double &k, double &b0,
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
    gsl_vector_set (x, 0, b0);

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
    fval = s->fval;

    b0 = gsl_vector_get (s->x, 0);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free(ss);
    return true;
}


template<typename TDOUBLE>
bool callGSL_simplex2(double &fval, double &tp, double &k, double &b0,
                  String<String<String<TDOUBLE> > > &statePosteriors, 
                  String<String<Observations> > &setObs, 
                  double &kMin, double &kMax,
                  AppOptions const& options)
{
    if (options.verbosity >= 2) 
        std::cout << "Call GSL multimin solver nmsimplex2 ..." << std::endl;

    int status;
    int iter = 0;
    int max_iter = options.maxIter_simplex;
    const size_t n = 2; 
    double size;
    double penalty = 0.01;  // fraction of function value*(-1) 

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    
    struct Params3<TDOUBLE> params = {tp, kMin, kMax, penalty, statePosteriors, setObs, options};
    gsl_multimin_function f;

    // instantiation of functor with all fixed params
    Fct_GSL_X_GAMMA2<TDOUBLE> fct(tp, kMin, kMax, penalty, statePosteriors, setObs, options);

    /* Set initial step sizes to 0.0001 */
    gsl_vector *ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 0.001);  // TODO different for k and theta (make sure not getting below 0!)
    // TODO adjust to given value 

    f.n = n;
    f.f = &fct_GSL_X_GAMMA2_W<TDOUBLE>;        // pointer to wrapper member function
    f.params =  &fct;       // pointer to functor (instead of to params)
    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, k);
    gsl_vector_set (x, 1, b0);
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
    fval = s->fval;

    // if k < kMin: fix k and optimize only for b0 (mean=exp(b0))
    if (gsl_vector_get (s->x, 0) < kMin)
    {
        std::cout << "Note: fixed shape parameter k to: " << kMin <<  std::endl;

        b0 = gsl_vector_get (s->x, 1);
        callGSL_simplex2_fixK(status, fval, tp, kMin, b0, statePosteriors, setObs, options);  

        gsl_vector_set (s->x, 0, kMin);
        gsl_vector_set (s->x, 1, b0);
    }
    else if (gsl_vector_get (s->x, 0) > kMax)
    {
        std::cout << "Note: fixed shape parameter k to: " << kMax << std::endl;
        b0 = gsl_vector_get (s->x, 1);
        callGSL_simplex2_fixK(status, fval, tp, kMax, b0, statePosteriors, setObs, options);  

        gsl_vector_set (s->x, 0, kMax);
        gsl_vector_set (s->x, 1, b0);
    }

    if (options.verbosity >= 2)
    {
        printf ("status = %s\n", gsl_strerror (status));
        std::cout << "GSL simplex2 .... k = " << gsl_vector_get (s->x, 0)  << " b0 = " << gsl_vector_get (s->x, 1) << std::endl;
    }
    k = gsl_vector_get (s->x, 0);
    b0 = gsl_vector_get (s->x, 1); 

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
    double fval = 100000.0;  // note: f was negated before, we minimze
    return callGSL_simplex2(fval, this->tp, this->k, this->b0, statePosteriors, setObs, kMin, kMax, options);    
}

template<typename TDOUBLE>
bool GAMMA2<TDOUBLE>::updateThetaAndK(String<String<double> > &startSet,
                    String<String<String<TDOUBLE> > > &statePosteriors, 
                    String<String<Observations> > &setObs, 
                    double &kMin, double &kMax,
                    AppOptions const&options)
{
    std::cout << "updateThetaAndK... kMax: " << kMax << std::endl;

    String<double> fvals;
    String<double> ks;
    String<double> b0s;
    resize(fvals, length(startSet), 100000.0, Exact()); // note: f was negated before, we minimze
    resize(ks, length(startSet), Exact());
    resize(b0s, length(startSet), Exact());

    // use multidimensional minimzation
    for (unsigned i = 0; i < length(startSet); ++i)
    {
        double k = startSet[i][0];  
        double b0 = startSet[i][1];
        if(!callGSL_simplex2(fvals[i], this->tp, k, b0, statePosteriors, setObs, kMin, kMax, options))
        {
            std::cout << "ERROR: during simplex optimization!" << std::endl;
            return false;
        }
        ks[i] = k;
        b0s[i] = b0;
    }
    double min_fval = 100000.0;
    for (unsigned i = 0; i < length(startSet); ++i)
    {
        if (fvals[i] < min_fval)
        {
            min_fval = fvals[i];
            this->k = ks[i];
            this->b0 = b0s[i];
        }
    }
    return true;  
}

//
///////////////////////////////////////////////


template<typename TDOUBLE>
long double GAMMA2<TDOUBLE>::getDensity(double const &x)   
{
    if (x < this->tp) return 0.0; 
    long double theta = (long double)exp(this->b0)/(long double)this->k;

    long double f1 = pow((long double)x, (long double)this->k - 1.0) * exp(-(long double)x/(long double)theta);
    long double f2 = pow((long double)theta, (long double)this->k) * tgamma((long double)this->k);
    
    // normalized lower incomplete gamma function
    long double nligf = boost::math::gamma_p((long double)this->k, (long double)this->tp/(long double)theta);
    if (nligf == 1.0) 
    {
        std::cout << "ERROR: (1 - nligf) is 0! Not set to max. value, should not happen for non-GLM gamma model!" << std::endl;
    } 
    return  ((f1/f2)/(1.0 - nligf));
}



//////////////////////////
// utils


template<typename TDOUBLE>
void myPrint(GAMMA2<TDOUBLE> &gamma)
{
    std::cout << "*** GAMMA2 ***" << std::endl;
    std::cout << "    b0:"<< gamma.b0 << std::endl;
    std::cout << "    mean:"<< exp(gamma.b0) << std::endl;
    std::cout << "    k:" << gamma.k << std::endl;
    std::cout << "    tp:" << gamma.tp << std::endl;
    std::cout << std::endl;
}

template<typename TDOUBLE>
bool checkConvergence(GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2, AppOptions &options)
{
    if (std::fabs(gamma1.b0 - gamma2.b0) > options.gamma_b_conv) return false;
    if (std::fabs(gamma1.k - gamma2.k) > options.gamma_k_conv) return false;

    return true;
}

template<typename TOut, typename TDOUBLE>
void printParams(TOut &out, GAMMA2<TDOUBLE> &gamma, int i)
{
    out << "gamma" << i << ".b0" << '\t' << gamma.b0 << std::endl;
    out << "gamma" << i << ".mean" << '\t' << exp(gamma.b0) << std::endl;
    out << "gamma" << i << ".k" << '\t' << gamma.k << std::endl;
    out << "gamma" << i << ".tp" << '\t' << gamma.tp << std::endl;
}


template<typename TDOUBLE>
void checkOrderG1G2(GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2, AppOptions &options)
{
    if (gamma1.k*gamma1.b0 > gamma2.b0)
    {
        std::swap(gamma1.b0, gamma2.b0);
        std::cout << "NOTE: swapped gamma1.b0 and gamma2.b0 ! " << std::endl;
    }
}


#endif
