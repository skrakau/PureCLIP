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

#ifndef APPS_HMMS_DENSITY_FUNCTIONS_CROSSLINK_REG_H_
#define APPS_HMMS_DENSITY_FUNCTIONS_CROSSLINK_REG_H_
   
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
#include <gsl/gsl_multiroots.h>

using namespace seqan;


////////
// ZTBIN_REG
////////
// P = (k-1)/(n-1) ?

class ZTBIN_REG
{
public:
    ZTBIN_REG(double b0_): b0(b0_) {}
    ZTBIN_REG() {}
 
    template<typename TType1, typename TType2> long double getDensity(TType1 const &k, TType2 const &n, double const &pred);
    template<typename TType1, typename TType2> long double getDensity(TType1 const &k, TType2 const &n);

    void updateP(String<String<String<long double> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const& options);
    void updateRegCoeffs(String<String<String<long double> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const&options);

    double b0;   // intercept
    String<double> regCoeffs;
};


// Functor for Brent's algorithm: find regression coefficients
// for given motif m; optimize b_m
struct FctLL_ZTBIN_REG
{
    FctLL_ZTBIN_REG(double const& b0_, char const& m_, String<String<String<long double> > > const& statePosteriors_,  String<String<Observations> > &setObs_, AppOptions const&options_) : b0(b0_), m(m_), statePosteriors(statePosteriors_), setObs(setObs_), options(options_)
    { 
    }
    double operator()(double const& b)
    {
        double ll = 0.0;       
        for (unsigned s = 0; s < 2; ++s)
        {
            String<double> lls;
            resize(lls, length(setObs[s]), 0.0, Exact());
#if HMM_PARALLEL
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)  
                {
                    // only optimize for sites with fimo score > 0.0 and corresponding to current id m
                    if (setObs[s][i].nEstimates[t] >= options.nThresholdForP && setObs[s][i].truncCounts[t] > 0 && setObs[s][i].motifIds[t] == m && setObs[s][i].fimoScores[t] > 0.0)
                    {
                        unsigned k = setObs[s][i].truncCounts[t];
                        unsigned n = (setObs[s][i].nEstimates[t] > setObs[s][i].truncCounts[t]) ? (setObs[s][i].nEstimates[t]) : (setObs[s][i].truncCounts[t]); 
                        double x = setObs[s][i].fimoScores[t];

                        if (((double)(k) / (double)(n)) <= options.maxkNratio)
                        {
                            double p = 1.0/(1.0+exp(-b0 - b*x));
                        
                            // l = log(1.0) -log(1.0 - pow((1.0-p), n)) + log (n over k) + k*log(p) + (n-k)*log(1.0-p);
                            // ignore parts not meaning for optimization! 
                            double l = -log(1.0 - pow((1.0-p), n)) + k*log(p) + (n-k)*log(1.0-p);
                            lls[i] += l * statePosteriors[s][i][t];
                        }
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
    double b0;
    char m;     // motif ID
    String<String<String<long double> > > statePosteriors;
    String<String<Observations> > &setObs;
    AppOptions options;
};


void ZTBIN_REG::updateRegCoeffs(String<String<String<long double> > > &statePosteriors, 
                         String<String<Observations> > &setObs, 
                         AppOptions const&options)
{ 
    int bits = 60;
    boost::uintmax_t maxIter = options.maxIter_brent;
    
    double bMin = 0.0;
    double bMax = 1.0;

    // for each input motif learn independent regCoeff (each position only one motif match with score assigned)
    for (unsigned char m = 0; m < options.nInputMotifs; ++m)
    {
        FctLL_ZTBIN_REG fct_ZTBIN_REG(this->b0, m, statePosteriors, setObs, options);
        std::pair<double, double> res = boost::math::tools::brent_find_minima(fct_ZTBIN_REG, bMin, bMax, bits, maxIter);         
        this->regCoeffs[m] = res.first;
    }
}



// use truncCounts
void ZTBIN_REG::updateP(String<String<String<long double> > > &statePosteriors, 
                  String<String<Observations> > &setObs, AppOptions const& options)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(setObs[s]); ++i)
        {
            for (unsigned t = 0; t < setObs[s][i].length(); ++t)
            {
                if (setObs[s][i].nEstimates[t] >= options.nThresholdForP && setObs[s][i].truncCounts[t] > 0 && setObs[s][i].fimoScores[t] == 0.0)      // avoid deviding by 0 (NOTE !), zero-truncated
                {
                    // p^ = (k-1)/(n-1); 'Truncated Binomial and Negative Binomial Distributions' Rider, 1955
                    unsigned k = setObs[s][i].truncCounts[t];
                    unsigned n = (setObs[s][i].nEstimates[t] > setObs[s][i].truncCounts[t]) ? (setObs[s][i].nEstimates[t]) : (setObs[s][i].truncCounts[t]);      
                    if (((double)(k) / (double)(n)) <= options.maxkNratio)
                    {
                        sum1 += statePosteriors[s][i][t] * ((double)(k - 1) / (double)(n - 1));        
                        sum2 += statePosteriors[s][i][t];
                    }
                 }
            }
        }
    }
    //std::cout << "updateP: sum1" << sum1 << " sum2: " << sum2 << " p: " << (sum1/sum2) << std::endl;
    double p = sum1/sum2;
    this->b0 = log(p/(1.0-p));

    updateRegCoeffs(statePosteriors, setObs, options);
}


// k: diagnostic events (de); n: read counts (c)
template<typename TType1, typename TType2> 
long double ZTBIN_REG::getDensity(TType1 const &k, TType2 const &n, double const &pred)
{
    if (k == 0) return 0.0;     // zero-truncated

    unsigned n2 = (n > k) ? (n) : (k);           // make sure n >= k      (or limit k?) 
   
    
    // use boost implementation, maybe avoids overflow
    boost::math::binomial_distribution<long double> boostBin;
    boostBin = boost::math::binomial_distribution<long double> ((int)n2, pred); 

    double res = boost::math::pdf(boostBin, k);
    if (std::isnan(res))   // or any other error?
    {
        std::cerr << "ERROR: binomial pdf is : " << res << std::endl;
        return 0.0;
    }
    return res * (1.0/(1.0 - pow((1.0 - pred), n2)));     // zero-truncated      TODO ???
}

template<typename TType1, typename TType2> 
long double ZTBIN_REG::getDensity(TType1 const &k, TType2 const &n)
{
    if (k == 0) return 0.0;     // zero-truncated

    unsigned n2 = (n > k) ? (n) : (k);          // make sure n >= k      (or limit k?) 
   
    long double pred = 1.0/(1.0+exp(- this->b0));

    // use boost implementation, maybe avoids overflow
    boost::math::binomial_distribution<long double> boostBin;
    boostBin = boost::math::binomial_distribution<long double> ((int)n2, pred); 

    long double res = boost::math::pdf(boostBin, k);
    if (std::isnan(res))   // or any other error?
    {
        std::cerr << "ERROR: binomial pdf is : " << res << std::endl;
        return 0.0;
    }
    return res * (1.0/(1.0 - pow((1.0 - pred), n2)));     // zero-truncated      TODO ???
}


void myPrint(ZTBIN_REG &bin)
{
    std::cout << "*** ZTBIN_REG ***" << std::endl;
    std::cout << "    p (0):"<< (1.0/(1.0+exp(-bin.b0))) << std::endl;
    std::cout << "    b0:"<< bin.b0 << std::endl;

    for (unsigned m = 0; m < length(bin.regCoeffs); ++m)
         std::cout << "    b" << m << ": " << bin.regCoeffs[m] << std::endl;

    std::cout << std::endl;
}


bool checkConvergence(ZTBIN_REG &bin1, ZTBIN_REG &bin2, AppOptions &options)
{
    if (std::fabs(bin1.b0 - bin2.b0) > options.bin_p_conv) return false;
    
    for (unsigned m = 0; m < length(bin1.regCoeffs); ++m)
        if (std::fabs(bin1.regCoeffs[m] - bin2.regCoeffs[m]) > options.bin_b_conv) return false;

    return true;
}


void checkOrderBin1Bin2(ZTBIN_REG &bin1, ZTBIN_REG &bin2)
{
    if (bin1.b0 > bin2.b0)
        std::swap(bin1.b0, bin2.b0); 
}


#endif
