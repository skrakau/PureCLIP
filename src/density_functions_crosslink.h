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


#ifndef APPS_HMMS_DENSITY_FUNCTIONS_CROSSLINK_H_
#define APPS_HMMS_DENSITY_FUNCTIONS_CROSSLINK_H_
   
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
// ZTBIN
////////
// P = (k-1)/(n-1) ?

class ZTBIN
{
public:
    ZTBIN(double p_): p(p_) {}
    ZTBIN() {}
 
    template<typename TType1, typename TType2> long double getDensity(TType1 const &k, TType2 const &n);

    void updateP(String<String<String<long double> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const& options); 

    double p;
};


// use truncCounts
void ZTBIN::updateP(String<String<String<long double> > > &statePosteriors, 
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
                if (setObs[s][i].nEstimates[t] >= options.nThresholdForP && setObs[s][i].truncCounts[t] > 0 )     // avoid deviding by 0 (NOTE !), zero-truncated
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
    this->p = sum1 / sum2;
}


// k: diagnostic events (de); n: read counts (c)
template<typename TType1, typename TType2> 
long double ZTBIN::getDensity(TType1 const &k, TType2 const &n)
{
    if (k == 0) return 0.0;     // zero-truncated

    unsigned n2 = (n > k) ? (n) : (k);          // make sure n >= k      (or limit k?) 
   
    // use boost implementation, maybe avoids overflow
    boost::math::binomial_distribution<long double> boostBin;
    boostBin = boost::math::binomial_distribution<long double> ((int)n2, this->p); 

    long double res = boost::math::pdf(boostBin, k);
    if (std::isnan(res))   // or any other error?
    {
        std::cerr << "ERROR: binomial pdf is : " << res << std::endl;
        return 0.0;
    }
    return res * (1.0/(1.0 - pow((1.0 - this->p), n2)));     // zero-truncated      TODO ???
}


void myPrint(ZTBIN &bin)
{
    std::cout << "*** ZTBIN ***" << std::endl;
    std::cout << "    p:"<< bin.p << std::endl;
    std::cout << std::endl;
}


bool checkConvergence(ZTBIN &bin1, ZTBIN &bin2, AppOptions &options)
{
    if (std::fabs(bin1.p - bin2.p) > options.bin_p_conv) return false;
    return true;
}


void checkOrderBin1Bin2(ZTBIN &bin1, ZTBIN &bin2)
{
    if (bin1.p > bin2.p)
        std::swap(bin1.p, bin2.p); 
}



#endif
