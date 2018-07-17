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

template<typename TDOUBLE>
class ZTBIN_REG
{
public:
    ZTBIN_REG(double b0_): b0(b0_) {}
    ZTBIN_REG() {}
 
    long double getDensity(unsigned const &k, unsigned const &n, double const &pred);
    long double getDensity(unsigned const &k, unsigned const &n);

    void updateP(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const& options);
    void updateRegCoeffs(String<String<String<TDOUBLE> > > &statePosteriors, String<String<Observations> > &setObs, AppOptions const&options);

    double b0;   // intercept
    String<double> regCoeffs;
};


// Functor for Brent's algorithm: find regression coefficients
// for given motif m; optimize b_m
template<typename TDOUBLE>
struct FctLL_ZTBIN_REG
{
    FctLL_ZTBIN_REG(double const& b0_, char const& m_, String<String<String<TDOUBLE> > > const& statePosteriors_,  String<String<Observations> > &setObs_, AppOptions const&options_) : b0(b0_), m(m_), statePosteriors(statePosteriors_), setObs(setObs_), options(options_)
    { 
    }
    double operator()(double const& b)
    {
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
                    // only optimize for sites with fimo score > 0.0 and corresponding to current id m
                    if (setObs[s][i].nEstimates[t] >= options.nThresholdForP && setObs[s][i].truncCounts[t] > 0 && setObs[s][i].motifIds[t] == m && setObs[s][i].fimoScores[t] > 0.0)
                    {
                        unsigned k = setObs[s][i].truncCounts[t];
                        unsigned n = (setObs[s][i].nEstimates[t] > setObs[s][i].truncCounts[t]) ? (setObs[s][i].nEstimates[t]) : (setObs[s][i].truncCounts[t]); 
                        double x = setObs[s][i].fimoScores[t];

                        if (((double)(k) / (double)(n)) <= options.maxkNratio)
                        {
                            TDOUBLE p = 1.0/(1.0+exp(-b0 - b*x));
                        
                            // l = log(1.0) -log(1.0 - pow((1.0-p), n)) + log (n over k) + k*log(p) + (n-k)*log(1.0-p);
                            // ignore parts not meaning for optimization! 
                            TDOUBLE l = -log(1.0 - pow((1.0-p), n)) + k*log(p) + (n-k)*log(1.0-p);
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
    String<String<String<TDOUBLE> > > statePosteriors;
    String<String<Observations> > &setObs;
    AppOptions options;
};


template<typename TDOUBLE>
void ZTBIN_REG<TDOUBLE>::updateRegCoeffs(String<String<String<TDOUBLE> > > &statePosteriors, 
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
        FctLL_ZTBIN_REG<TDOUBLE> fct_ZTBIN_REG(this->b0, m, statePosteriors, setObs, options);
        std::pair<double, double> res = boost::math::tools::brent_find_minima(fct_ZTBIN_REG, bMin, bMax, bits, maxIter);         
        this->regCoeffs[m] = res.first;
    }
}



// use truncCounts
template<typename TDOUBLE>
void ZTBIN_REG<TDOUBLE>::updateP(String<String<String<TDOUBLE> > > &statePosteriors, 
                  String<String<Observations> > &setObs, AppOptions const& options)
{
    TDOUBLE sum1 = 0.0;
    TDOUBLE sum2 = 0.0;
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
                        sum1 += statePosteriors[s][i][t] * ((TDOUBLE)(k - 1) / (TDOUBLE)(n - 1));        
                        sum2 += statePosteriors[s][i][t];
                    }
                 }
            }
        }
    }
    //std::cout << "updateP: sum1" << sum1 << " sum2: " << sum2 << " p: " << (sum1/sum2) << std::endl;
    TDOUBLE p = sum1/sum2;
    this->b0 = log(p/(1.0-p));

    updateRegCoeffs(statePosteriors, setObs, options);
}


// k: diagnostic events (de); n: read counts (c)
template<typename TDOUBLE> 
long double ZTBIN_REG<TDOUBLE>::getDensity(unsigned const &k, unsigned const &n, double const &pred)
{
    if (k == 0) return 0.0;     // zero-truncated

    unsigned n2 = (n > k) ? (n) : (k);           // make sure n >= k      (or limit k?) 
   
    
    // use boost implementation, maybe avoids overflow
    boost::math::binomial_distribution<TDOUBLE> boostBin;
    boostBin = boost::math::binomial_distribution<TDOUBLE> ((int)n2, pred); 

    TDOUBLE res = boost::math::pdf(boostBin, k);
    if (std::isnan(res))   // or any other error?
    {
        std::cerr << "ERROR: binomial pdf is : " << res << std::endl;
        return 0.0;
    }
    return res * (TDOUBLE)(1.0/(1.0 - pow((1.0 - pred), n2)));     // zero-truncated      TODO ???
}

template<typename TDOUBLE> 
long double ZTBIN_REG<TDOUBLE>::getDensity(unsigned const &k, unsigned const &n)
{
    if (k == 0) return 0.0;     // zero-truncated

    unsigned n2 = (n > k) ? (n) : (k);          // make sure n >= k      (or limit k?) 
   
    double pred = 1.0/(1.0+exp(- this->b0));

    // use boost implementation, maybe avoids overflow
    boost::math::binomial_distribution<TDOUBLE> boostBin;
    boostBin = boost::math::binomial_distribution<TDOUBLE> ((int)n2, pred); 

    TDOUBLE res = boost::math::pdf(boostBin, k);
    if (std::isnan(res))   // or any other error?
    {
        std::cerr << "ERROR: binomial pdf is : " << res << std::endl;
        return 0.0;
    }
    return res * (TDOUBLE)(1.0/(1.0 - pow((1.0 - pred), n2)));     // zero-truncated      TODO ???
}


//////////////////////////
// utils


template<typename TDOUBLE>
bool loadBinParams(ZTBIN_REG<TDOUBLE> &bin1, ZTBIN_REG<TDOUBLE> &bin2, AppOptions &options)
{
    std::ifstream ifs(toCString(options.inParFileName));
    if (!ifs.good())
        std::cerr << "ERROR: Could not open file containing input parameters!\n";
    std::string lineBuffer;

    // 1) get number of motif regression coefficients
    unsigned m1 = 0;
    unsigned m2 = 0;
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
                std::cerr << "ERROR: could not read first value\n";
                return false;
            }

            CharString str = value1.c_str();
            Prefix<CharString >::Type pref = prefix(str, 6);
            CharString suf = suffix(str, 6);

            int sufInt = std::atoi(toCString(suf));

            if (pref == "bin1.b" && sufInt > m1)
                m1 = sufInt;
            else if (pref == "bin2.b" && sufInt > m2)
                m2 = sufInt; 
        }
    }
    if (m1 != m2) 
    {
        std::cerr << "ERROR: not same number of regression coefficients given for bin1 and bin2!\n";
        return false;
    }
    resize(bin1.regCoeffs, m1, 0.0, Exact());
    resize(bin2.regCoeffs, m1, 0.0, Exact());

    // 2) go to beginning of file again and read line by line
    ifs.clear();
    ifs.seekg(0, std::ios::beg);
    while (std::getline(ifs, lineBuffer))
    {
        //std::cout << lineBuffer << std::endl;
        std::string value1;
        std::string value2;
        std::stringstream ss(lineBuffer);
        if (!ss.str().empty())
        {
            if (!std::getline(ss, value1, '\t'))
                std::cerr << "ERROR: could not read first value\n";

            if (value1.c_str() == "bin1.b0")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for bin1.b0\n";
                    return false;
                }
                bin1.b0 = std::strtod(value2.c_str(), NULL);
            }
            else if (value1.c_str() == "bin2.b0")
            {
                if (!std::getline(ss, value2, '\t')) 
                {
                    std::cerr << "ERROR: could not read second value for bin2.b0\n";
                    return false;
                }
                bin2.b0 = std::strtod(value2.c_str(), NULL);
            }   
            else
            {
                for (unsigned m = 0; m <= m1; ++m)
                {
                    std::stringstream test1; 
                    std::stringstream test2; 
                    test1 << "bin1.b" << m;
                    test2 << "bin2.b" << m;
                    if (value1 == test1.str())    //?
                    {
                        if (!std::getline(ss, value2, '\t')) 
                        {
                            std::cerr << "ERROR: could not read second value for bin1.b" << m << "\n";
                            return false;
                        }
                        bin1.regCoeffs[m] = std::strtod(value2.c_str(), NULL);
                        break;
                    }
                    else if (value1 == test2.str())   // ?
                    {
                        if (!std::getline(ss, value2, '\t')) 
                        {
                            std::cerr << "ERROR: could not read second value for bin2.b" << m << "\n";
                            return false;
                        }
                        bin1.regCoeffs[m] = std::strtod(value2.c_str(), NULL);
                        break;
                    }  
                }
            }
        }
    }
    return true;
}


template<typename TDOUBLE>
void myPrint(ZTBIN_REG<TDOUBLE> &bin)
{
    std::cout << "*** ZTBIN_REG ***" << std::endl;
    std::cout << "    p (0):"<< (1.0/(1.0+exp(-bin.b0))) << std::endl;
    std::cout << "    b0:"<< bin.b0 << std::endl;

    for (unsigned m = 0; m < length(bin.regCoeffs); ++m)
         std::cout << "    b" << (m+1) << ": " << bin.regCoeffs[m] << std::endl;

    std::cout << std::endl;
}


template<typename TOut, typename TDOUBLE>
void printParams(TOut &out, ZTBIN_REG<TDOUBLE> &bin, int i)
{
    out << "bin" << i << ".b0" << '\t' << bin.b0 << std::endl;

    for (unsigned m = 0; m < length(bin.regCoeffs); ++m)
        out << "bin" << i << ".b" << (m+1) << '\t' << bin.regCoeffs[m] << std::endl;
    out << std::endl;    
}


template<typename TDOUBLE>
bool checkConvergence(ZTBIN_REG<TDOUBLE> &bin1, ZTBIN_REG<TDOUBLE> &bin2, AppOptions &options)
{
    if (std::fabs(bin1.b0 - bin2.b0) > options.bin_p_conv) return false;
    
    for (unsigned m = 0; m < length(bin1.regCoeffs); ++m)
        if (std::fabs(bin1.regCoeffs[m] - bin2.regCoeffs[m]) > options.bin_b_conv) return false;

    return true;
}


template<typename TDOUBLE>
void checkOrderBin1Bin2(ZTBIN_REG<TDOUBLE> &bin1, ZTBIN_REG<TDOUBLE> &bin2)
{
    if (bin1.b0 > bin2.b0)
        std::swap(bin1.b0, bin2.b0); 
}


#endif
