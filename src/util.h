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


#ifndef APPS_HMMS_UTIL_H_
#define APPS_HMMS_UTIL_H_

#include <iostream>
#include <fstream>
#include <seqan/bed_io.h>

#include <math.h>    

using namespace seqan;

namespace seqan {

    // TODO investigate impact
    class LogSumExp_lookupTable
    {
        public:
            unsigned size;
            double minValue;
            String<long double> lookupTable;

            LogSumExp_lookupTable() : size(0), minValue(0.0) {}

            LogSumExp_lookupTable(unsigned size_, double minValue_) : size(size_), minValue(minValue_) 
        {
            resize(lookupTable, size+1, Exact());
            for(int i = 0; i <= size; ++i)
                lookupTable[i] = log1p(exp(i*(-minValue/size) + minValue));
            std::cout << "Created look-up table for values from " << ((0 * -minValue/size) + minValue) << " to " << ((size * -minValue/size) + minValue) << " with step size " << (-minValue/size) << " (size: " << size<< ")." << std::endl;
        }

            long double logSumExp_add(long double f1, long double f2) const
            {
                if (f1 > f2)
                {
                    if (f2 - f1 < minValue)
                    {
                        //std::cout << "Note: f1 - f2 below precision! " << (f2 - f1) << " (returning larger value)." <<  std::endl;
                        return f1;
                    }
                    return f1 + lookupTable[(int)(((f2-f1) - minValue) * (size/-minValue))];
                }
                else
                {
                    if (f1 - f2 < minValue)
                    {
                        //std::cout << "Note: f2 - f1 below precision! " << (f1 - f2) << " (returning larger value)." <<  std::endl;
                        return f2;
                    }
                    return f2 + lookupTable[(int)(((f1-f2) - minValue) * (size/-minValue))];
                }
            }
    };    

    struct AppOptions
    {
        CharString bamFileName;
        CharString baiFileName;
        CharString refFileName;
        CharString outFileName;
        CharString outRegionsFileName;
        CharString parFileName;
        CharString rpkmFileName;
        CharString inputBamFileName;
        CharString inputBaiFileName;
        CharString fimoFileName;

	    bool crosslinkAtTruncSite;
        CharString                  intervals_str;
        String<unsigned>            intervals_contigIds;
        String<String<unsigned> >   intervals_positions;
        CharString                  applyChr_str;
        String<unsigned>            applyChr_contigIds;

        bool posteriorDecoding;
        double prior_kdeThreshold;
        unsigned prior_enrichmentThreshold;
        unsigned maxIter_brent;
        unsigned maxIter_bw;
        unsigned maxIter_simplex;
        double g1_kMin;
        double g2_kMin;
        double g1_kMax;
        double g2_kMax;
        bool g1_k_le_g2_k;
        double p1;
        double p2;
        double gamma_k_conv;
        double gamma_theta_conv;
        double gamma_b_conv;
        double bin_p_conv;
        double bin_b_conv;
        unsigned binSize;
        unsigned bandwidth;
        unsigned bandwidthN;
        unsigned nKernelGap;
        unsigned intervalOffset;

        bool gaussianKernel;
        bool epanechnikovKernel;
        double useKdeThreshold;

        bool estimateNfromKdes;
        double nThresholdForP;
        double nThresholdForTransP;
        double get_nThreshold;
        double minTransProbCS;
        double maxkNratio;
        unsigned maxBinN;

        unsigned polyAThreshold;
        bool excludePolyAFromLearning;
        bool excludePolyTFromLearning;
        bool excludePolyA;
        bool excludePolyT;

        bool gslSimplex2;
        long double min_nligf;
        double kMin_simplex;
        double kMax_simplex;

        bool useCov_RPKM;
        bool useLogRPKM;
        double minRPKMtoFit;
        bool mrtf_kdeSglt;
        bool discardSingletonIntervals;
        unsigned maxTruncCount;
        unsigned maxTruncCount2;

        bool useFimoScore;
        unsigned nInputMotifs;

        unsigned distMerge;
        bool useHighPrecision;      // long double to compute emission probabilities, state posteriors, Forward-Backward (alpha, beta) values
        LogSumExp_lookupTable lookUp;   // table containing log-sum-exp precomputed values to avoid expensive log and exp operations
        unsigned lookupTable_size;
        double lookupTable_minValue;
        unsigned selectRead;

        unsigned numThreads;
        unsigned numThreadsA;
        bool outputAll;
        // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
        int verbosity;

        AppOptions() :
            crosslinkAtTruncSite(false),
            posteriorDecoding(true),
            prior_enrichmentThreshold(7),   // KDE threshold is used corresponding to 7 read starts at one position
            maxIter_brent(100),              // brent
            maxIter_bw(50),                  // baum-welch
            maxIter_simplex(2000),           // simplex            
            g1_kMin(1.0),                   // shape parameter for gamma distribution; set min. to avoid eProbs getting zero!
            g2_kMin(1.0),
            g1_kMax(10.0),
            g2_kMax(10.0),
            g1_k_le_g2_k(true),
            p1(0.01),                       // initial values for bin1.p
            p2(0.15),                       // .. bin2.p
            gamma_k_conv(0.0001),
            gamma_theta_conv(0.0001),
            gamma_b_conv(0.0001),
            bin_p_conv(0.0001),
            bin_b_conv(0.0001),
            binSize(0),                     // if not specified: 2* bdw
            bandwidth(50),                  // h, standard deviation for gaussian kernel
            bandwidthN(0),                  // .... used for estionation of N
            nKernelGap(0),                  // 
            intervalOffset(50),             // offset for covered intervals to be stored in observations
            gaussianKernel(true),
            epanechnikovKernel(false),
            useKdeThreshold(0.0),
            estimateNfromKdes(true),
            nThresholdForP(10),             // threshold regarding n used for fitting p, if GLM, this need to be larger!
            nThresholdForTransP(0),         // threshold regarding n used for fitting p, if GLM, this need to be larger! (usually state post. for 'enriched' enough) 
            get_nThreshold(false),          // estimate threshold based on expected read start counts
            minTransProbCS(0.0001),
            maxkNratio(1.0),                // ignore sites for binomial learning with ratio greater (maybe caused by mapping artifacts)
            maxBinN(50000),                 // sites above not used for learning
            polyAThreshold(10),
            excludePolyAFromLearning(false),
            excludePolyTFromLearning(false),
            excludePolyA(false),
            excludePolyT(false),
            gslSimplex2(true),
            min_nligf(0.99999999),          // min. normalized lower incomplete gamma function. NOTE: precission of boost computation is limited, set to min. value in order to avoid 1s!
            kMin_simplex(0.5),              // not used currently ...
            kMax_simplex(15.0),
            useCov_RPKM(false),
            useLogRPKM(true),
            minRPKMtoFit(-5.0),
            mrtf_kdeSglt(true),                 // use singleton KDE value as mrtf for GLM fitting (assuming same bandwidth for input KDEs!)
            discardSingletonIntervals(true),    // delete intervals with singleton reads to save memory (and runtime) !! influence on transProbs?
            maxTruncCount(500),                 // used to ignore intervals for learning
            maxTruncCount2(65000),              // to store, larger values are truncated to this value, avoid overflow of n 
            useFimoScore(false),
            nInputMotifs(1),
            distMerge(8),
            useHighPrecision(false),
            lookupTable_size(600000),
            lookupTable_minValue(-2000.0),
            selectRead(0),
            numThreads(1),
            numThreadsA(0),
            outputAll(false),
            verbosity(1)
        {}
    };

    class Times
    {
    public:
        double time_all;
        double time_loadObservations;
        double time_learnHMM;
        double time_applyHMM;
        double time_applyHMM2;      // in parallel

        static Times & instance()
        {
            static Times times;
            return times;
        }

    private:
        Times() :
            time_all(0),
            time_loadObservations(0),
            time_learnHMM(0),
            time_applyHMM(0),
            time_applyHMM2(0)
        {}
    };


    struct ContigObservations {
        String<__uint16> truncCounts; 
    };

    void reverse(ContigObservations &contigObservations)
    {    
        reverse(contigObservations.truncCounts); 
    }


    // workaround because partially specialized member function are forbidden
    // wrapper class for observations
    struct Observations {
        Infix<String<__uint16> >::Type truncCounts;
        unsigned contigId; 

        String<__uint32>    nEstimates;      
        String<double>      kdes;       // used for 'enriched'. 'non-enriched' classification
        String<double>      kdesN;    // used to estimate the binomial n parameters (decoupled, might be useful e.g. for longer crosslink clusters) 
        String<double>      rpkms;      // TODO change name -> e.g. bgSignal
        String<float>       fimoScores; // for each t: one motif score
        String<char>        motifIds; // for each t: one motif score
        bool                discard;    // NOTE: only use for application, not for learning!

        Observations(Infix<String<__uint16> >::Type _truncCounts) : truncCounts(_truncCounts),
                                                                    discard(false) {}
        Observations() : truncCounts(),
                         discard(false) {}

        void estimateNs(AppOptions &options);                       // using raw counts
        void estimateNs(double b0, double b1, AppOptions /*&options*/); // using KDEs
        void computeKDEs(AppOptions &options);
        void computeKDEs(String<__uint16> &inputTruncCounts, AppOptions &options);    // input signal

        unsigned length();
    };

    inline unsigned Observations::length() {
        return endPosition(this->truncCounts) - beginPosition(this->truncCounts); // length(this->truncCounts);
    }

    
    void Observations::estimateNs(AppOptions &options)  
    { 
        resize(this->nEstimates, length(), Exact());
        unsigned w_50 = floor((double)options.bandwidthN - 0.1);    // binSize should be odd

        for (unsigned t = 0; t < length(); ++t)
        {
            unsigned sum = 0;
            for (unsigned j = std::max((int)t - (int)w_50, (int)0); (j < length()) && (j <= t + w_50); ++j)  // inefficient... update on the fly
                sum += this->truncCounts[j];
            
             this->nEstimates[t] = std::max(sum, (unsigned)1); 
        }
    }

    // use simple linear regression, estimate from KDE values
    void Observations::estimateNs(double b0, double b1, AppOptions /*&options*/)    
    { 
        resize(this->nEstimates, length(), Exact());
        for (unsigned t = 0; t < length(); ++t)
        {
            this->nEstimates[t] = (__uint16)std::max((int)floor(b0 + b1*this->kdesN[t] + 0.5), 1);   // avoid becoming 0  
        }
    }


    ///////////////////////////////
    // gaussian kernel for smothing
    template<typename TType> 
    double getGaussianKernelDensity(TType const &u)
    {
        double fac1 = 1.0 / (double)std::sqrt(2.0*M_PI);
        double fac2 = (-1.0)*pow(u, 2);

        return (fac1 * exp(fac2/2.0));
    }
    ///////////////////////////////
    //  Epanechnikov kernel for smothing
    template<typename TType> 
    double getEpanechnikovKernelDensity(TType const &u)
    {
        if (std::abs(u) > 1) return 0.0;
        
        return (3.0/4.0 * (1.0 - pow(u, 2)));
    }
    ///////////////////////////////
    // custom kernel for estimation of N
    // lower weight for very close by neighbour positions (other crosslinks within same motif regions for example)
    template<typename TType> 
    double getCustomKernelDensity(TType const &u, AppOptions &options)
    {
        TType cu = u;
        // simply cutoff
        if (u*options.bandwidthN <= options.nKernelGap && u*options.bandwidthN != 0) return cu = (double)options.nKernelGap/(double)options.bandwidthN;       

        double fac1 = 1.0 / (double)std::sqrt(2.0*M_PI);
        double fac2 = (-1.0)*pow(cu, 2);


        if (u*options.bandwidthN == 0) return (fac1 * exp(fac2/2.0));     
        if (u*options.bandwidthN <= options.nKernelGap) return (fac1 * exp(fac2/2.0))*0.0;     

        return (fac1 * exp(fac2/2.0));
    }

    void Observations::computeKDEs(AppOptions &options)
    {
        resize(this->kdes, length());

        unsigned w_50 = options.bandwidth * 4;
        // precompute kernel densities
        String<double> kernelDensities;
        resize(kernelDensities, w_50 + 1, 0.0);

        // precompute kernel densities   -> K(d/h) store at position d
        for (unsigned i = 0; i <= w_50; ++i)
        {
            if (options.gaussianKernel)
                kernelDensities[i] = getGaussianKernelDensity((double)i/(double)options.bandwidth);
            else if (options.epanechnikovKernel)
                kernelDensities[i] = getEpanechnikovKernelDensity((double)i/(double)options.bandwidth);        
        }
        for (unsigned t = 0; t < length(); ++t)
        {
            double kde = 0.0;
            for (unsigned i = std::max((int)t - (int)w_50, (int)0); (i < length()) && (i <= t + w_50); ++i)  // inefficient if low genome coverage and not selected only for covered regions!!!
            { 
                kde += this->truncCounts[i] * kernelDensities[(unsigned)std::abs((int)t - (int)i)];
            }
            this->kdes[t] = kde/(double)options.bandwidth; 
        }

        /////////////////////////////////
        // same, but for kdes used for estimation of n
        // TODO problem: interval size dependent on main bandwidth parameter, 
        /////////////////////////////////
        resize(this->kdesN, length());

        w_50 = options.bandwidthN * 4;
        // precompute kernel densities
        clear(kernelDensities);
        resize(kernelDensities, w_50 + 1, 0.0);

        // precompute kernel densities   -> K(d/h) store at position d
        for (unsigned i = 0; i <= w_50; ++i)
        {
            kernelDensities[i] = getCustomKernelDensity((double)i/(double)options.bandwidthN, options);      
        }
        for (unsigned t = 0; t < length(); ++t)
        {
            double kde = 0.0;
            for (unsigned i = std::max((int)t - (int)w_50, (int)0); (i < length()) && (i <= t + w_50); ++i)  // inefficient if low genome coverage and not selected only for covered regions!!!
            { 
                kde += this->truncCounts[i] * kernelDensities[(unsigned)std::abs((int)t - (int)i)];
            }
            this->kdesN[t] = kde/(double)options.bandwidthN; 
        }
    }

    // for input truncCounts (not stored in observations)
    // keep in mind: same intervals used as for target (+- 2*bdw)
    // could cause underestimation of input KDEs at interval boarders
    // anyway only very low values of gaussian kernel there
    void Observations::computeKDEs(String<__uint16> &truncCounts, AppOptions &options)
    {
        resize(this->rpkms, length());

        unsigned w_50 = options.bandwidth * 4;
        // precompute kernel densities
        String<double> kernelDensities;
        resize(kernelDensities, w_50 + 1, 0.0);

        // precompute kernel densities   -> K(d/h) store at position d
        for (unsigned i = 0; i <= w_50; ++i)
        {
            if (options.gaussianKernel)
                kernelDensities[i] = getGaussianKernelDensity((double)i/(double)options.bandwidth);
            else if (options.epanechnikovKernel)
                kernelDensities[i] = getEpanechnikovKernelDensity((double)i/(double)options.bandwidth);        
        }
        for (unsigned t = 0; t < length(); ++t)
        {
            double kde = 0.0;
            for (unsigned i = std::max((int)t - (int)w_50, (int)0); (i < length()) && (i <= t + w_50); ++i)  // inefficient if low genome coverage and not selected only for covered regions!!!
            { 
                kde += truncCounts[i] * kernelDensities[(unsigned)std::abs((int)t - (int)i)];
            }
            if (options.useLogRPKM)
            {
                if ((kde/(double)options.bandwidth) > 0.0)
                    this->rpkms[t] = log(kde/(double)options.bandwidth); 
                else 
                    this->rpkms[t] = options.minRPKMtoFit - 1.0; 
            }
            else
                this->rpkms[t] = kde/(double)options.bandwidth; 
        }
    }


    struct Data {
        String<String<Observations> >               setObs;       // F/R:interval:t
        String<String<unsigned> >                   setPos;
        String<String<String<String<double> > > >   statePosteriors;  // F/R:state:interval:t
        String<String<String<__uint8> > >           states;
    };

    void append(Data &dataA, Data &dataB)
    {
        for (unsigned s = 0; s < 2; ++s)
        {
            if (!empty(dataB.setObs[s]))
                append(dataA.setObs[s], dataB.setObs[s]);
            if (!empty(dataB.setPos[s]))        
                append(dataA.setPos[s], dataB.setPos[s]);
            if (!empty(dataB.statePosteriors[s]))
                append(dataA.statePosteriors[s], dataB.statePosteriors[s]);
            if (!empty(dataB.states[s]))
                append(dataA.states[s], dataB.states[s]); 
        }
    }

    void clear(Data &data)
    {
        clear(data.setObs);
        clear(data.setPos);
        clear(data.statePosteriors);
        clear(data.states);
    }


}




#endif
