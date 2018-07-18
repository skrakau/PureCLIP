//    ======================================================================
//    PureCLIP: capturing target-specific protein-RNA interaction footprints
//    ======================================================================
//    Copyright (C) 2017  Sabrina Krakau, Max Planck Institute for Molecular 
//    Genetics
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef APPS_HMMS_DENSITY_FUNCTIONS_REG_H_
#define APPS_HMMS_DENSITY_FUNCTIONS_REG_H_
   
#include <iostream>
#include <fstream>
#include <math.h>       // lgamma 

#include <boost/math/tools/minima.hpp>      // BRENT's algorithm
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>       // normalized lower incomplete gamma function: gamma_p()
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
//#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>

using namespace seqan;



//////////////////////////////////////////////////////////////////////////
// GAMMA2_REG: left threshold, forced to be zero
//////////////////////////////////////////////////////////////////////////

template<typename TDOUBLE>
class GAMMA2_REG  // ignore positions with KDE below theshold
{
public:
    GAMMA2_REG(double mean_, double k_, double tp_): mean(mean_), k(k_), tp(tp_) {}
    GAMMA2_REG(double tp_): tp(tp_) {}
    GAMMA2_REG() {}

    long double getDensity(double const &kde, double const &pred);
    void updateMean(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const& options); 
    void updateK(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, double &kMin, double &kMax, AppOptions const& options);

    bool updateRegCoeffsAndK(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, double &kMin, double &kMax, AppOptions const& options); 
 

    double mean;   
    double b0;
    double b1;
    double k;       // shape parameter 
    double tp;      // truncation point
};



////////////////////////////////////
// GSL newton 

template<typename TDOUBLE>
struct Fct_GSL_N_GAMMA2_REG
{
    Fct_GSL_N_GAMMA2_REG(double const & tp_, double const& k_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_,  
                                  String<String<Observations> > & setObs_,  
                                  AppOptions const&options_) : tp(tp_), k(k_),
                                                               statePosteriors(statePosteriors_),  
                                                               setObs(setObs_), 
                                                               options(options_)
    { 
    }
    // f
    double operator()(const gsl_vector * x)
    {       
        const double b0 = gsl_vector_get (x, 0);
        const double b1 = gsl_vector_get (x, 1);

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
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1 && setObs[s][i].rpkms[t] >= options.minRPKMtoFit)
                    {
                        double kde = setObs[s][i].kdes[t];
                        double x1 = setObs[s][i].rpkms[t];
                        double pred = exp(b0 + b1 * x1);

                        double nligf = boost::math::gamma_p(k, (tp*k/pred));
            
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

    // f'
    void Gradient(const gsl_vector * x,  gsl_vector * g) {

        const double b0 = gsl_vector_get (x, 0);
        const double b1 = gsl_vector_get (x, 1);

        TDOUBLE f_0 = 0.0;
        TDOUBLE f_1 = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<TDOUBLE> f_0_S;
            String<TDOUBLE> f_1_S;
            resize(f_0_S, length(setObs[s]), 0.0, Exact());
            resize(f_1_S, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)
                {    
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1 && setObs[s][i].rpkms[t] >= options.minRPKMtoFit)
                    {
                        double kde = setObs[s][i].kdes[t];
                        double x1 = setObs[s][i].rpkms[t];
                        double pred  = exp(b0 + b1 * x1);

                        double h = tp*k/pred;
                        double ligf = boost::math::tgamma_lower(k, h);

                        // b0
                        TDOUBLE p_0 = k*kde/pred - k - exp(-h)*pow(h, k)/(tgamma(k) - ligf);   
                        // b1
                        TDOUBLE p_1 = x1*p_0;   

                        f_0_S[i] +=  p_0 * statePosteriors[s][i][t];
                        f_1_S[i] +=  p_1 * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                f_0 += f_0_S[i];
                f_1 += f_1_S[i];
            }
        }

        gsl_vector_set (g, 0, -f_0);
        gsl_vector_set (g, 1, -f_1);
    }
    // f and f''
    void FdF(const gsl_vector * x, double * f, gsl_vector * g) {

        const double b0 = gsl_vector_get (x, 0);
        const double b1 = gsl_vector_get (x, 1);

        *f = 0.0;
        TDOUBLE f_0 = 0.0;
        TDOUBLE f_1 = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<TDOUBLE> f_S;
            String<TDOUBLE> f_0_S;
            String<TDOUBLE> f_1_S;
            resize(f_S, length(setObs[s]), 0.0, Exact());
            resize(f_0_S, length(setObs[s]), 0.0, Exact());
            resize(f_1_S, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)
                {    
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1 && setObs[s][i].rpkms[t] >= options.minRPKMtoFit)
                    {
                        double kde = setObs[s][i].kdes[t];
                        double x1 = setObs[s][i].rpkms[t];
                        double pred = exp(b0 + b1 * x1);

                        double h = tp*k/pred;
                        double nligf = boost::math::gamma_p(k, h);
                        double ligf = boost::math::tgamma_lower(k, h);
            
                        TDOUBLE p = (k-1.0)*log(kde) - k * (kde/pred + log(pred)) - k*log(1.0/k) - lgamma(k) - log(1.0 - nligf);
                        f_S[i] +=  p * statePosteriors[s][i][t];

                        // b0
                        TDOUBLE p_0 = k*kde/pred - k - exp(-h)*pow(h, k)/(tgamma(k) - ligf);   
                        // b1
                        TDOUBLE p_1 = x1*p_0;   
                        f_0_S[i] +=  p_0 * statePosteriors[s][i][t];
                        f_1_S[i] +=  p_1 * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                *f += f_S[i];
                f_0 += f_0_S[i];
                f_1 += f_1_S[i];
            }
        }
        *f *= (-1);
        gsl_vector_set (g, 0, -f_0);
        gsl_vector_set (g, 1, -f_1);
   }

private:
    double tp;
    double k;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > & setObs;
    AppOptions options;
};


// Wrapper functions for functors
template<typename TDOUBLE>
double fct_GSL_N_GAMMA2_REG_f_W (const gsl_vector * x, void * p) {

    Fct_GSL_N_GAMMA2_REG<TDOUBLE> * function = reinterpret_cast< Fct_GSL_N_GAMMA2_REG<TDOUBLE> *> (p);
    return (*function)( x );        
} 

template<typename TDOUBLE>
void fct_GSL_N_GAMMA2_REG_df_W (const gsl_vector * x, void * p,  gsl_vector * g) {

    Fct_GSL_N_GAMMA2_REG<TDOUBLE> * function = reinterpret_cast< Fct_GSL_N_GAMMA2_REG<TDOUBLE> *> (p);
    (*function).Gradient( x, g );
}

template<typename TDOUBLE>
void fct_GSL_N_GAMMA2_REG_fdf_W (const gsl_vector * x, void * p, double *f, gsl_vector * g ) {

    Fct_GSL_N_GAMMA2_REG<TDOUBLE> * function = reinterpret_cast< Fct_GSL_N_GAMMA2_REG<TDOUBLE> *> (p);
    (*function).FdF( x, f, g);
} 


template<typename TDOUBLE>
struct Params2
{
    double tp;
    double k;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > setObs;
    AppOptions options;
};



void print_state2(size_t iter, gsl_multimin_fdfminimizer * s, void * /*fct*/)
{
    gsl_vector *gradients = gsl_multimin_fdfminimizer_gradient (s);
    
    printf ("iter = %5lu b0 = % 10.7f b1 = % 10.7f  "
            "g0 = % 10.7f g1 = % 10.7f "
            "f(x) = % 10.7f \n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (gradients, 0),
            gsl_vector_get (gradients, 1),
            gsl_multimin_fdfminimizer_minimum (s));

    // test if left and right higher than minimum
    /*double test_f;
    gsl_vector *test_g = gsl_vector_alloc (2);

    gsl_vector *test_x = gsl_vector_alloc (2);
    gsl_vector_set (test_x, 0, gsl_vector_get (s->x, 0));
    gsl_vector_set (test_x, 1, gsl_vector_get (s->x, 1));

    double new_x = gsl_vector_get (test_x, 0) + 0.01;
    gsl_vector_set (test_x, 0, new_x);

    fct_GSL_N_GAMMA2_REG_fdf_W(test_x, fct, &test_f, test_g);
    std::cout << "test: left: f(): " << test_f << '\t';

    new_x = gsl_vector_get (test_x, 0) - 0.02;
    gsl_vector_set (test_x, 0, new_x);

    fct_GSL_N_GAMMA2_REG_fdf_W(test_x, fct, &test_f, test_g);
    std::cout << "right f(): " << test_f << std::endl;

    gsl_vector_free (test_x);
    gsl_vector_free (test_g);
    */
}


template<typename TDOUBLE>
int callGSL_newton2(double &tp, double &k, double &b0, double &b1,
                  String<String<String<TDOUBLE> > > &statePosteriors, 
                  String<String<Observations> > &setObs, 
                  AppOptions const& options)
{
    std::cout << "Call GSL multiroot solver ..." << std::endl;
    int status;
    int iter = 0;
    int max_iter = options.maxIter_simplex;
    const size_t n = 2; 

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    struct Params2<TDOUBLE> params = {tp, k, statePosteriors, setObs, options};
    gsl_multimin_function_fdf f;

    // instantiation of functor with all fixed params
    Fct_GSL_N_GAMMA2_REG<TDOUBLE> fct(tp, k, statePosteriors, setObs, options);

    f.n = n;
    f.f = &fct_GSL_N_GAMMA2_REG_f_W<TDOUBLE>;        // pointer to wrapper member function
    f.df = &fct_GSL_N_GAMMA2_REG_df_W<TDOUBLE>;
    f.fdf = &fct_GSL_N_GAMMA2_REG_fdf_W<TDOUBLE>;
    f.params =  &fct;       // pointer to functor (instead of to params)

    gsl_vector *x = gsl_vector_alloc (n);

    gsl_vector_set (x, 0, b0);
    gsl_vector_set (x, 1, b1);

    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, n);
    gsl_multimin_fdfminimizer_set (s, &f, x, 0.0001, 1e-4);   // initial step size, line minimization parameter

    print_state2 (iter, s, &fct);

    do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);  

      if (status)
        break;

      gsl_vector *gradients = gsl_multimin_fdfminimizer_gradient (s);
      status = gsl_multimin_test_gradient (gradients, 1e-4);

      if (status == GSL_SUCCESS)
        printf ("Minimum found at:\n");

      print_state2 (iter, s, &fct);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    printf ("status = %s\n", gsl_strerror (status));

    std::cout << "GSL newton .... b0 = " << gsl_vector_get (s->x, 0)  << " b1 = " << gsl_vector_get (s->x, 1) << std::endl;
    b0 = gsl_vector_get (s->x, 0);
    b1 = gsl_vector_get (s->x, 1); 

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
    return 0;
}


template<typename TDOUBLE>
void GAMMA2_REG<TDOUBLE>::updateMean(String<String<String<TDOUBLE> > > &statePosteriors, 
                    String<String<Observations> > &setObs,  
                    AppOptions const&options)
{
    // use multidimensional rootfinding, newton
    callGSL_newton2(this->tp, this->k, this->b0, this->b1, statePosteriors, setObs, options);
}



//////////////////////////////////////////////
// update betas and k together using simplex2
//////////////////////////////////////////////

template<typename TDOUBLE>
struct Fct_GSL_X_GAMMA2_REG
{
    Fct_GSL_X_GAMMA2_REG(double const & tp_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_,
                                  String<String<Observations> > & setObs_,  
                                  AppOptions const&options_) : tp(tp_),
                                                               statePosteriors(statePosteriors_),  
                                                               setObs(setObs_), 
                                                               options(options_)
    { 
    }
    // f
    double operator()(const gsl_vector * x)
    {      
        const double k = gsl_vector_get (x, 0);
        const double b0 = gsl_vector_get (x, 1);
        const double b1 = gsl_vector_get (x, 2);

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
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1 && setObs[s][i].rpkms[t] >= options.minRPKMtoFit)
                    {
                        double kde = setObs[s][i].kdes[t];
                        double x1 = setObs[s][i].rpkms[t];
                        double pred = exp(b0 + b1 * x1);

                        double nligf = boost::math::gamma_p(k, (tp*k/pred));
            
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
struct Fct_GSL_X_GAMMA2_REG_fixK
{
    Fct_GSL_X_GAMMA2_REG_fixK(double const & tp_, double const & k_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_,
                                  String<String<Observations> > & setObs_,
                                  AppOptions const&options_) : tp(tp_), k(k_),
                                                               statePosteriors(statePosteriors_),
                                                               setObs(setObs_),  
                                                               options(options_)
    { 
    }
    // f
    double operator()(const gsl_vector * x)
    {      
        const double b0 = gsl_vector_get (x, 0);
        const double b1 = gsl_vector_get (x, 1);

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
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1 && setObs[s][i].rpkms[t] >= options.minRPKMtoFit)
                    {
                        double kde = setObs[s][i].kdes[t];
                        double x1 = setObs[s][i].rpkms[t];
                        double pred = exp(b0 + b1 * x1);

                        double nligf = boost::math::gamma_p(k, (tp*k/pred));
            
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
    String<String<Observations> > & setObs;
    AppOptions options;
};


// Wrapper functions for functors
template<typename TDOUBLE>
double fct_GSL_X_GAMMA2_REG_W (const gsl_vector * x, void * p) {

    Fct_GSL_X_GAMMA2_REG<TDOUBLE> * function = reinterpret_cast< Fct_GSL_X_GAMMA2_REG<TDOUBLE> *> (p);
    return (*function)( x );        
}

template<typename TDOUBLE>
double fct_GSL_X_GAMMA2_REG_fixK_W (const gsl_vector * x, void * p) {

    Fct_GSL_X_GAMMA2_REG_fixK<TDOUBLE> * function = reinterpret_cast< Fct_GSL_X_GAMMA2_REG_fixK<TDOUBLE> *> (p);
    return (*function)( x );        
} 

template<typename TDOUBLE>
struct Params5
{
    double tp;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > setObs;
    AppOptions options;
};

template<typename TDOUBLE>
bool callGSL_simplex2_fixK(int &status, 
                  double &tp, double &k, double &b0, double &b1,
                  String<String<String<TDOUBLE> > > &statePosteriors, 
                  String<String<Observations> > &setObs, 
                  AppOptions const& options)
{
    int iter = 0;
    int max_iter = options.maxIter_simplex;
    const size_t n = 2; 
    double size;

    std::cout << "Note: Set k to: " << k << std::endl;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    
    struct Params2<TDOUBLE> params = {tp, k, statePosteriors, setObs, options};
    gsl_multimin_function f;

    // instantiation of functor with all fixed params
    Fct_GSL_X_GAMMA2_REG_fixK<TDOUBLE> fct(tp, k, statePosteriors, setObs, options);

    /* Set initial step sizes to */
    gsl_vector *ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 0.001);  

    f.n = n;
    f.f = &fct_GSL_X_GAMMA2_REG_fixK_W<TDOUBLE>;        // pointer to wrapper member function
    f.params =  &fct;       // pointer to functor (instead of to params)

    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, b0);
    gsl_vector_set (x, 1, b1);

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

            printf ("%5d %10.7f %10.7f f() = %7.7f size = %.7f\n", 
                  iter,
                  gsl_vector_get (s->x, 0), 
                  gsl_vector_get (s->x, 1), 
                  s->fval, size);
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    b0 = gsl_vector_get (s->x, 0);
    b1 = gsl_vector_get (s->x, 1);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free(ss);
    return true;
}

template<typename TDOUBLE>
bool callGSL_simplex2(double &tp, double &k, double &b0, double &b1,
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
    const size_t n = 3; 
    double size;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    
    struct Params5<TDOUBLE> params = {tp, statePosteriors, setObs, options};
    gsl_multimin_function f;

    // instantiation of functor with all fixed params
    Fct_GSL_X_GAMMA2_REG<TDOUBLE> fct(tp, statePosteriors, setObs, options);

    /* Set initial step sizes to 0.0001 */
    gsl_vector *ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 0.001);  // TODO different for k and theta (make sure not getting below 0!)
    // TODO adjust to given value 

    f.n = n;
    f.f = &fct_GSL_X_GAMMA2_REG_W<TDOUBLE>;        // pointer to wrapper member function
    f.params =  &fct;       // pointer to functor (instead of to params)

    gsl_vector *x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, k);
    gsl_vector_set (x, 1, b0);
    gsl_vector_set (x, 2, b1);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, n);
    gsl_multimin_fminimizer_set (s, &f, x, ss);   

    if (options.verbosity >= 2)
    {
        printf ("%5d %10.7f %10.7f %10.7f", 
                  0,
                  gsl_vector_get (s->x, 0), 
                  gsl_vector_get (s->x, 1),
                  gsl_vector_get (s->x, 2));
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
        if (gsl_vector_get (s->x, 0) < kMin)
        {
            std::cout << "Note: limited shape parameter k to: " << kMin << std::endl;   //". This could be caused by outliers: high peaks, potentially background binding. Check if transcripts/chromosomes used for learning are representative." <<  std::endl;

            b0 = gsl_vector_get (s->x, 1);
            b1 = gsl_vector_get (s->x, 2);
            callGSL_simplex2_fixK(status, tp, kMin, b0, b1, statePosteriors, setObs, options);  

            gsl_vector_set (s->x, 0, kMin);
            gsl_vector_set (s->x, 1, b0);
            gsl_vector_set (s->x, 2, b1);
            break;
        }
        else if (gsl_vector_get (s->x, 0) > kMax)
        {
            std::cout << "Note: limited shape parameter k to: " << kMax << std::endl; 

            b0 = gsl_vector_get (s->x, 1);
            b1 = gsl_vector_get (s->x, 2);
            callGSL_simplex2_fixK(status, tp, kMax, b0, b1, statePosteriors, setObs, options);  

            gsl_vector_set (s->x, 0, kMax);
            gsl_vector_set (s->x, 1, b0);
            gsl_vector_set (s->x, 2, b1);
           break;
        }

        if (options.verbosity >= 2)
        {
            printf ("%5d %10.7f %10.7f %10.7f f() = %7.7f size = %.7f\n", 
                  iter,
                  gsl_vector_get (s->x, 0), 
                  gsl_vector_get (s->x, 1),
                  gsl_vector_get (s->x, 2),
                  s->fval, size);
        }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    if (options.verbosity >= 2)
    {
        printf ("status = %s\n", gsl_strerror (status));
        std::cout << "GSL simplex2 .... k = " << gsl_vector_get (s->x, 0)  << " b0 = " << gsl_vector_get (s->x, 1) << " b1 = " << gsl_vector_get (s->x, 2) << std::endl;
    }

    k = gsl_vector_get (s->x, 0);
    b0 = gsl_vector_get (s->x, 1); 
    b1 = gsl_vector_get (s->x, 2); 

    if (b1 < 0.0) 
    {
        std::cerr << "ERROR: b1 became < 0! Should be >= 0." << std::endl;
        return false;
    }

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free(ss);
    return true;
}



template<typename TDOUBLE>
bool GAMMA2_REG<TDOUBLE>::updateRegCoeffsAndK(String<String<String<TDOUBLE> > > &statePosteriors, 
                    String<String<Observations> > &setObs,  
                    double &kMin, double &kMax,
                    AppOptions const&options)
{
    // use multidimensional rootfinding, newton
    return callGSL_simplex2(this->tp, this->k, this->b0, this->b1, statePosteriors, setObs, kMin, kMax, options);
}




///////////////////////////////////////
// update k
///////////////////////////////////////


/*void GAMMA2_REG::approximateK(String<String<double> > &statePosteriorsF, 
                              String<String<double> > &statePosteriorsR, 
                              String<Observations> &setObsF, String<Observations> &setObsR, 
                              AppOptions const&options)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    // forward
    for (unsigned i = 0; i < length(setObsF); ++i)
    {
        for (unsigned t = 0; t < setObsF[i].length(); ++t)
        {
            if (setObsF[i].kdes[t] >= options.useKdeThreshold && setObsF[i].kdes[t] >= this->tp && setObsF[i].truncCounts[t] >= 1 && setObsF[i].rpkms[t] >= options.minRPKMtoFit)
            {
                sum1 += statePosteriorsF[i][t] * setObsF[i].kdes[t];
                sum2 += statePosteriorsF[i][t] * log(setObsF[i].kdes[t]);
                sum3 += statePosteriorsF[i][t];
            }
        }
    }
    // reverse
    for (unsigned i = 0; i < length(setObsR); ++i)
    {
        for (unsigned t = 0; t < setObsR[i].length(); ++t)
        {
            if (setObsR[i].kdes[t] >= options.useKdeThreshold && setObsR[i].kdes[t] >= this->tp && setObsR[i].truncCounts[t] >= 1 && setObsR[i].rpkms[t] >= options.minRPKMtoFit)
            {
                sum1 += statePosteriorsR[i][t] * setObsR[i].kdes[t];
                sum2 += statePosteriorsR[i][t] * log(setObsR[i].kdes[t]);
                sum3 += statePosteriorsR[i][t];
            }
        }
    }
    double s = log(sum1/sum3) - (sum2/sum3);

    this->k = (3.0 - s + sqrt(pow(s - 3.0, 2) + 24.0*s))/(12.0 * s);
}*/



// Functor for Brent's algorithm
// maximize for k
template<typename TDOUBLE>
struct Fct_GAMMA2_REG_k
{
    Fct_GAMMA2_REG_k(double const& b0_, double const& b1_, 
                                  String<String<String<TDOUBLE> > > const& statePosteriors_, 
                                  String<String<Observations> > & setObs_, 
                                  AppOptions const&options_) : b0(b0_), b1(b1_), 
                                                               statePosteriors(statePosteriors_),  
                                                               setObs(setObs_), 
                                                               options(options_)
    { 
    }
    double operator()(double const& k)
    {
        // Group log-likelihood function evaluations regarding binned kde vaues ! todo

        TDOUBLE ll = 0.0;
        for (unsigned s = 0; s < 2; ++s)
        {
            String<TDOUBLE> llsS;
            resize(llsS, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif 
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)  
                {
                    if (setObs[s][i].kdes[t] >= options.useKdeThreshold && setObs[s][i].truncCounts[t] >= 1 && 
                        setObs[s][i].rpkms[t] >= options.minRPKMtoFit)
                    {
                        double kde = setObs[s][i].kdes[t];
                        double x1 = setObs[s][i].rpkms[t];
                        double pred = exp(b0 + b1 * x1);

                        double nligf = boost::math::gamma_p(k, (options.useKdeThreshold/(pred/k)));

                        TDOUBLE p = (k-1.0)*log(kde) - k * (kde/pred + log(pred)) - k*log(1.0/k) - lgamma(k);
                        p -= log(1.0 - nligf);
                        llsS[i] +=  p * statePosteriors[s][i][t];
                    }
                }
            }
            // combine results from threads
            for (unsigned i = 0; i < length(setObs[s]); ++i)
                ll += llsS[i];
        }
        return (-ll);
    }

private:
    double b0;
    double b1;
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > & setObs;
    AppOptions options;
};


template<typename TDOUBLE>
void GAMMA2_REG<TDOUBLE>::updateK(String<String<String<TDOUBLE> > > &statePosteriors, 
                         String<String<Observations> > &setObs, 
                         double &kMin, double &kMax,
                         AppOptions const&options)
{ 
    int bits = 60;
    boost::uintmax_t maxIter = options.maxIter_brent;
    
    Fct_GAMMA2_REG_k<TDOUBLE> fct_GAMMA2_REG_k(this->b0, this->b1, statePosteriors, setObs, options);
    std::pair<double, double> res = boost::math::tools::brent_find_minima(fct_GAMMA2_REG_k, kMin, kMax, bits, maxIter);         // use somehow initial guess to save time? or interval around prev. value?

    this->k = res.first;
}


template<typename TDOUBLE>
long double GAMMA2_REG<TDOUBLE>::getDensity(double const &kde, double const &pred)   
{
    if (kde < this->tp) return 0.0;

    //double pred = exp(this->b0 + this->b1 * x);

    double theta = pred/this->k;
    // if (kde == 0.0) should not occur, checked while computing eProbs
    TDOUBLE f1 = pow(kde, this->k - 1.0) * exp(-kde/theta);
    TDOUBLE f2 = pow(theta, this->k) * tgamma(this->k);
    if (f2 ==  0.0) std::cout << "ERROR: f2 is 0!" << std::endl;


    // normalized lower incomplete gamma function
    double nligf = boost::math::gamma_p(this->k, this->tp/theta);
    if ((1.0 - nligf) == 0.0) std::cout << "ERROR: (1 - nligf) is 0!"  << " kde: " << kde << " pred: " << pred << std::endl;

    return  ((f1/f2)/(1.0 - nligf));
}






//////////////////////////
// utils


template<typename TDOUBLE>
void myPrint(GAMMA2_REG<TDOUBLE> &gamma)
{
    std::cout << "*** GAMMA2_REG ***" << std::endl;
    std::cout << "    b0:"<< gamma.b0 << std::endl;
    std::cout << "    b1:"<< gamma.b1 << std::endl;
    std::cout << "    mean:"<< gamma.mean << std::endl;
    std::cout << "    k:" << gamma.k << std::endl;
    std::cout << "    tp:" << gamma.tp << std::endl;
    std::cout << std::endl;
}

template<typename TDOUBLE>
bool checkConvergence(GAMMA2_REG<TDOUBLE> &gamma1, GAMMA2_REG<TDOUBLE> &gamma2, AppOptions &options)
{
    if (std::fabs(gamma1.b0 - gamma2.b0) > options.gamma_b_conv) return false;
    if (std::fabs(gamma1.b1 - gamma2.b1) > options.gamma_b_conv) return false;
    if (std::fabs(gamma1.k - gamma2.k) > options.gamma_k_conv) return false;

    return true;
}

template<typename TOut, typename TDOUBLE>
void printParams(TOut &out, GAMMA2_REG<TDOUBLE> &gamma, int i)
{
    out << "gamma" << i << ".b0" << '\t' << gamma.b0 << std::endl;
    out << "gamma" << i << ".b1" << '\t' << gamma.b1 << std::endl;
    out << "gamma" << i << ".theta" << '\t' << exp(gamma.b0)/gamma.k << std::endl;
    out << "gamma" << i << ".k" << '\t' << gamma.k << std::endl;
    out << "gamma" << i << ".tp" << '\t' << gamma.tp << std::endl;
    out << std::endl;    
}


template<typename TDOUBLE>
void checkOrderG1G2(GAMMA2_REG<TDOUBLE> &gamma1, GAMMA2_REG<TDOUBLE> &gamma2, AppOptions &options)
{
    if (exp(gamma1.b0) > exp(gamma2.b0))
    {
        std::swap(gamma1.b0, gamma2.b0);
        std::cout << "NOTE: swapped gamma1.b0 and gamma2.b0 ! " << std::endl;
    }
    if (options.g1_k_le_g2_k && gamma1.k > gamma2.k)
    {
        gamma2.k = gamma1.k;
        std::cout << "NOTE: set gamma2.k to gamma1.k ! " << std::endl;
    }

}


#endif
