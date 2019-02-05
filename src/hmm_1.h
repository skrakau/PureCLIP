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

#ifndef APPS_HMMS_HMM_1_H_
#define APPS_HMMS_HMM_1_H_
  

#include <iostream>
#include <fstream>

#include "density_functions.h"
#include "density_functions_reg.h"
#include "density_functions_crosslink.h"
#include "density_functions_crosslink_reg.h"
#include <math.h>  

using namespace seqan;


template <typename TGAMMA, typename TBIN>
class HMM {     

public:

    __uint8                                 K;                  // no. of sates
    String<String<String<long double> > >   initProbs;          // initial probabilities

    String<String<Observations> >           & setObs;          // workaround for partial specialization
    String<String<unsigned> >               & setPos;
    unsigned                                contigLength;
    String<String<long double> >            transMatrix;

    HMM(int K_, String<String<Observations> > & setObs_, String<String<unsigned> > & setPos_, unsigned &contigLength_): K(K_), setObs(setObs_), setPos(setPos_), contigLength(contigLength_)
    {
        // initialize transition probabilities
        resize(transMatrix, K, Exact());
        for (unsigned i = 0; i < K; ++i)
            resize(transMatrix[i], K, 0.25, Exact());
       
        resize(initProbs, 2, Exact());
        resize(eProbs, 2, Exact());
        resize(statePosteriors, 2, Exact());
        for (unsigned s = 0; s < 2; ++s)
        {
            resize(initProbs[s], length(setObs[s]), Exact());
            resize(eProbs[s], length(setObs[s]), Exact());
            resize(statePosteriors[s], K, Exact());
            for (unsigned k = 0; k < K; ++k)
                resize(statePosteriors[s][k], length(setObs[s]), Exact());

            for (unsigned i = 0; i < length(setObs[s]); ++i)
            {
                // set initial probabilities to uniform
                resize(initProbs[s][i], K, Exact());
                for (unsigned k = 0; k < K; ++k)
                    initProbs[s][i][k] = 1.0/K;

                unsigned T = setObs[s][i].length();
                resize(eProbs[s][i], T, Exact());
                for (unsigned k = 0; k < K; ++k)
                    resize(statePosteriors[s][k][i], T, Exact());

                for (unsigned t = 0; t < T; ++t)
                {
                    resize(eProbs[s][i][t], K, Exact());
                }
            }
        }
    }
    // Copy constructor 
//     HMM<TGAMMA, TBIN>(const HMM<TGAMMA, TBIN> &hmm2) {x = p2.x; y = p2.y; } 

    HMM<TGAMMA, TBIN>();
    ~HMM<TGAMMA, TBIN>();
    
    bool computeEmissionProbs(ModelParams<TGAMMA, TBIN> &modelParams, bool learning, AppOptions &options);
    bool iForward(String<String<long double> > &alphas_1, unsigned s, unsigned i, String<String<long double> > &logA, AppOptions &options);    
    bool iBackward(String<String<long double> > &betas_1, unsigned s, unsigned i, String<String<long double> > &logA, AppOptions &options);    
    bool computeStatePosteriorsFB(AppOptions &options);
    bool computeStatePosteriorsFBupdateTrans(AppOptions &options);
    bool updateTransAndPostProbs(AppOptions &options);
    bool updateDensityParams(TGAMMA &gamma1, TGAMMA &gamma2, unsigned &iter, unsigned &trial, AppOptions &options);
    bool updateDensityParams(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &options);
    bool baumWelch(ModelParams<TGAMMA, TBIN> &modelParams, CharString learnTag, AppOptions &options);
    bool applyParameters(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &/*options*/);
    void posteriorDecoding(String<String<String<__uint8> > > &states);
    void rmBoarderArtifacts(String<String<String<__uint8> > > &states, String<Data> &data_replicates, String<ModelParams<TGAMMA, TBIN> > &modelParams);

    // for each F/R,interval,t, state ....
    String<String<String<String<double> > > > eProbs;           // emission/observation probabilities  P(Y_t | S_t) -> precompute for each t given Y_t = (C_t, T_t) !!!
    String<String<String<String<double> > > > statePosteriors;  // for each k: for each covered interval string of posteriors
};


template<typename TGAMMA, typename TBIN>
HMM<TGAMMA, TBIN>::~HMM<TGAMMA, TBIN>()
{
    clear(this->eProbs);
    clear(this->statePosteriors);
    clear(this->initProbs);
    clear(this->transMatrix);
   // do not touch observations
}


/////////////////////////////////////////////////////////////////
// functionalities for computations in log space
/////////////////////////////////////////////////////////////////

long double myLog(long double x)
{
    if (x == 0) return std::numeric_limits<long double>::quiet_NaN();  
    return log(x);
}

long double myExp(long double x)
{
    if (std::isnan(x)) return 0.0;  
    return exp(x);
}

// log-sum-exp trick
long double get_logSumExp(long double &f1, long double &f2, LogSumExp_lookupTable &lookUp)
{
    if (std::isnan(f1)) return f2;
    if (std::isnan(f2)) return f1;

    if (std::isinf(f1)) return f1;
    if (std::isinf(f2)) return f2;

    return lookUp.logSumExp_add(f1, f2);
}

// log-sum-exp trick
long double get_logSumExp_states(long double f1, long double f2, long double f3, long double f4, LogSumExp_lookupTable &lookUp)
{
    long double sum; 
    sum = get_logSumExp(f1, f2, lookUp);
    sum = get_logSumExp(sum, f3, lookUp);
    sum = get_logSumExp(sum, f4, lookUp);
    return sum;
}

// log-sum-exp trick for string
long double get_logSumExp(String<long double> &fs, LogSumExp_lookupTable &lookUp)
{   
    long double sum = std::numeric_limits<long double>::quiet_NaN();

    for (unsigned i = 0; i < length(fs); ++i)
        sum = get_logSumExp(sum, fs[i], lookUp);

    return sum;
}

// log-sum-exp trick for String of String
long double get_logSumExp(String<String<long double> > &fs, LogSumExp_lookupTable &lookUp)
{
    long double sum = std::numeric_limits<long double>::quiet_NaN();

    for (unsigned i = 0; i < length(fs); ++i)
        for (unsigned j = 0; j < length(fs[i]); ++j)
            sum = get_logSumExp(sum, fs[i][j], lookUp);

    return sum;
}


/////////////////////////////////////////////////////////////////
// emission probabilities
/////////////////////////////////////////////////////////////////


template<typename TEProbs, typename TSetObs>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA &gamma1, GAMMA &gamma2, ZTBIN &bin1, ZTBIN &bin2, unsigned t, AppOptions &options)
{
    long double gamma1_d = 1.0;
    long double gamma2_d = 0.0;
    if (setObs.kdes[t] >= gamma1.tp) 
    {
        gamma1_d = gamma1.getDensity(setObs.kdes[t]);
        gamma2_d = gamma2.getDensity(setObs.kdes[t]); 
    }
    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
    }
    // log-space
    eProbs[0] = myLog(gamma1_d) + myLog(bin1_d); 
    eProbs[1] = myLog(gamma1_d) + myLog(bin2_d);
    eProbs[2] = myLog(gamma2_d) + myLog(bin1_d);
    eProbs[3] = myLog(gamma2_d) + myLog(bin2_d);

    // check if valid
    if ((gamma1_d + gamma2_d == 0.0) || (bin1_d + bin2_d == 0.0) ||
       (std::isnan(eProbs[0]) && std::isnan(eProbs[1]) && std::isnan(eProbs[2]) && std::isnan(eProbs[3])) )
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "WARNING: emission probabilities 0.0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << gamma1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << gamma2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 0.0;
        eProbs[1] = std::numeric_limits<double>::quiet_NaN();
        eProbs[2] = std::numeric_limits<double>::quiet_NaN();
        eProbs[3] = std::numeric_limits<double>::quiet_NaN();
        return false;
    } 
    return true;
}

template<typename TEProbs, typename TSetObs>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA_REG &gamma1, GAMMA_REG &gamma2, ZTBIN &bin1, ZTBIN &bin2, unsigned t, AppOptions &options)
{
    long double x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
    long double gamma1_pred = exp(gamma1.b0 + gamma1.b1 * x);
    long double gamma2_pred = exp(gamma2.b0 + gamma2.b1 * x);

    long double gamma1_d = 1.0;
    long double gamma2_d = 0.0;
    if (setObs.kdes[t] >= gamma1.tp) 
    {
        gamma1_d = gamma1.getDensity(setObs.kdes[t], gamma1_pred, options);
        gamma2_d = gamma2.getDensity(setObs.kdes[t], gamma2_pred, options); 
    }

    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
    }
    // log-space
    eProbs[0] = myLog(gamma1_d) + myLog(bin1_d); 
    eProbs[1] = myLog(gamma1_d) + myLog(bin2_d);
    eProbs[2] = myLog(gamma2_d) + myLog(bin1_d);
    eProbs[3] = myLog(gamma2_d) + myLog(bin2_d);

    // 
    if ((gamma1_d + gamma2_d == 0.0) || (bin1_d + bin2_d == 0.0) ||
            (std::isnan(eProbs[0]) && std::isnan(eProbs[1]) && std::isnan(eProbs[2]) && std::isnan(eProbs[3])) )
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "WARNING: emission probabilities going against 0.0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate b: " << x << " predicted mean 'non-enriched': " << gamma1_pred << " predicted mean 'enriched': " << gamma2_pred << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << gamma1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << gamma2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }    
        eProbs[0] = 0.0;
        eProbs[1] = std::numeric_limits<double>::quiet_NaN();
        eProbs[2] = std::numeric_limits<double>::quiet_NaN();
        eProbs[3] = std::numeric_limits<double>::quiet_NaN();
        return false;
    }
    return true;
}

template<typename TEProbs, typename TSetObs>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA &gamma1, GAMMA &gamma2, ZTBIN_REG &bin1, ZTBIN_REG &bin2, unsigned t, AppOptions &options)
{
    long double gamma1_d = 1.0;
    long double gamma2_d = 0.0;
    if (setObs.kdes[t] >= gamma1.tp) 
    {
        gamma1_d = gamma1.getDensity(setObs.kdes[t]);
        gamma2_d = gamma2.getDensity(setObs.kdes[t]); 
    }
    unsigned mId = setObs.motifIds[t];
    long double bin1_pred = 1.0/(1.0+exp(-bin1.b0 - bin1.regCoeffs[mId]*setObs.fimoScores[t]));
    long double bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*setObs.fimoScores[t]));

    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred, options);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred, options);
    }

    eProbs[0] = myLog(gamma1_d) + myLog(bin1_d);
    eProbs[1] = myLog(gamma1_d) + myLog(bin2_d);
    eProbs[2] = myLog(gamma2_d) + myLog(bin1_d);
    eProbs[3] = myLog(gamma2_d) + myLog(bin2_d);

    //
    if ((gamma1_d + gamma2_d == 0.0) || (bin1_d + bin2_d == 0.0) ||
            (std::isnan(eProbs[0]) && std::isnan(eProbs[1]) && std::isnan(eProbs[2]) && std::isnan(eProbs[3])) )
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "WARNING: emission probabilities going against 0.0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate x: " << setObs.fimoScores[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << gamma1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << gamma2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 0.0;
        eProbs[1] = std::numeric_limits<double>::quiet_NaN();
        eProbs[2] = std::numeric_limits<double>::quiet_NaN();
        eProbs[3] = std::numeric_limits<double>::quiet_NaN();
        return false;
    }
    return true;
}

template<typename TEProbs, typename TSetObs>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA_REG &gamma1, GAMMA_REG &gamma2, ZTBIN_REG &bin1, ZTBIN_REG &bin2, unsigned t, AppOptions &options)
{
    long double x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
    long double gamma1_pred = exp(gamma1.b0 + gamma1.b1 * x);
    long double gamma2_pred = exp(gamma2.b0 + gamma2.b1 * x);

    long double gamma1_d = 1.0;
    long double gamma2_d = 0.0;
    if (setObs.kdes[t] >= gamma1.tp) 
    {
        gamma1_d = gamma1.getDensity(setObs.kdes[t], gamma1_pred, options);
        gamma2_d = gamma2.getDensity(setObs.kdes[t], gamma2_pred, options); 
    }
    unsigned mId = setObs.motifIds[t];
    long double bin1_pred = 1.0/(1.0+exp(-bin1.b0 - bin1.regCoeffs[mId]*setObs.fimoScores[t]));
    long double bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*setObs.fimoScores[t]));

    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred, options);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred, options);
    }

    eProbs[0] = myLog(gamma1_d) + myLog(bin1_d);
    eProbs[1] = myLog(gamma1_d) + myLog(bin2_d);
    eProbs[2] = myLog(gamma2_d) + myLog(bin1_d);
    eProbs[3] = myLog(gamma2_d) + myLog(bin2_d);

    // 
    if ((gamma1_d + gamma2_d == 0.0) || (bin1_d + bin2_d == 0.0) ||
            (std::isnan(eProbs[0]) && std::isnan(eProbs[1]) && std::isnan(eProbs[2]) && std::isnan(eProbs[3])) )
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "WARNING: emission probabilities going against 0.0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate b: " << x << " predicted mean 'non-enriched': " << gamma1_pred << " predicted mean 'enriched': " << gamma2_pred << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate x: " << setObs.fimoScores[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << gamma1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << gamma2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 0.0;
        eProbs[1] = std::numeric_limits<double>::quiet_NaN();
        eProbs[2] = std::numeric_limits<double>::quiet_NaN();
        eProbs[3] = std::numeric_limits<double>::quiet_NaN();
        return false;
    }
    return true;
}


template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::computeEmissionProbs(ModelParams<TGAMMA, TBIN> &modelParams, bool learning, AppOptions &options)
{
    bool stop = false;
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            bool discardInterval = false;
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)  
            {
                if (this->setObs[s][i].kdes[t] == 0.0)
                {
                    std::cerr << "ERROR: KDE is 0.0 at i " << i << " t: " << t << std::endl;
                    SEQAN_OMP_PRAGMA(critical) 
                    stop = true;
                }

                if (!computeEProb(this->eProbs[s][i][t], this->setObs[s][i], modelParams.gamma1, modelParams.gamma2, modelParams.bin1, modelParams.bin2, t, options))
                {
                    SEQAN_OMP_PRAGMA(critical) 
                    discardInterval = true;
                }
            }
            if (learning && discardInterval)
            {
                SEQAN_OMP_PRAGMA(critical) 
                std::cout << "ERROR: Emission probability became 0.0! This might be due to artifacts or outliers." << std::endl;
                SEQAN_OMP_PRAGMA(critical)
                if (options.verbosity >= 2)
                {
                    if (s == 0) 
                        std::cout << " Interval: [" << (this->setPos[s][i]) << ", " << (this->setPos[s][i] + this->setObs[s][i].length()) << ") on forward strand." << std::endl;
                    else 
                        std::cout << " Interval: [" << (this->contigLength - this->setPos[s][i] - 1) << ", " << (this->contigLength - this->setPos[s][i] - 1 + this->setObs[s][i].length()) << ") on reverse strand." << std::endl;
                }
                stop = true;
                if (!options.useHighPrecision)  // TODO ?
                {
                    SEQAN_OMP_PRAGMA(critical) 
                    std::cout << "NOTE: Try running PureCLIP in high floating-point precision mode (long double, parameter '-ld')." << std::endl;
                }
            }
            else if (!learning && discardInterval) 
            {
                this->setObs[s][i].discard = true;
                SEQAN_OMP_PRAGMA(critical) 
                std::cout << "Warning: discarding interval on forward strand due to emission probabilities of 0.0 (set to state 'non-enriched + non-crosslink')." << std::endl;
                if (options.verbosity >= 2)
                {
                    SEQAN_OMP_PRAGMA(critical) 
                    if (s == 0) 
                        std::cout << " Interval [" << (this->setPos[s][i]) << ", " << (this->setPos[s][i] + this->setObs[s][i].length()) << ") on forward strand. " << std::endl;
                    else 
                        std::cout << " Interval [" << (this->contigLength - this->setPos[s][i] - 1) << ", " << (this->contigLength - this->setPos[s][i] - 1 + this->setObs[s][i].length()) << ") on reverse strand." << std::endl;
                }
                if (!options.useHighPrecision)  // TODO ?
                {
                    SEQAN_OMP_PRAGMA(critical) 
                    std::cout << "NOTE: If this happens frequently, rerun PureCLIP in high floating-point precision mode (long double, parameter '-ld')." << std::endl;
                }
            }
        }   
    }
    if (stop) return false;
    return true;
}



/////////////////////////////////////////////////////////////////
// forward-backward algorithm parts
/////////////////////////////////////////////////////////////////

// for one interval only
// Forward-algorithm: log-space
template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::iForward(String<String<long double> > &alphas_1, unsigned s, unsigned i, String<String<long double> > &logA, AppOptions &options)
{
    // NOTE
    // in log-space: alphas_1, eProbs
    // trans. probs, init porbs, state post. probs. not in log-space

    // for t = 1
    for (unsigned k = 0; k < this->K; ++k)
    {
        alphas_1[0][k] = myLog(this->initProbs[s][i][k]) + this->eProbs[s][i][0][k];    // log ? initProbs should not become 0.0!
    }

    // for t = 2:T
    for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
    {
        for (unsigned k = 0; k < this->K; ++k)
        { 
            long double f1 = alphas_1[t-1][0] + logA[0][k] + this->eProbs[s][i][t][k];
            long double f2 = alphas_1[t-1][1] + logA[1][k] + this->eProbs[s][i][t][k];
            long double f3 = alphas_1[t-1][2] + logA[2][k] + this->eProbs[s][i][t][k];
            long double f4 = alphas_1[t-1][3] + logA[3][k] + this->eProbs[s][i][t][k];

            alphas_1[t][k] = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);

            if (std::isinf(alphas_1[t][k]))
            {
                std::cout << "ERROR: alphas_1[" << t << "][" << k << "] is " << alphas_1[t][k] << std::endl;
                std::cout << "       f1 " << f1 << " f2 " << f2 << " f3 " << f3 << " f4 " << f4 << std::endl;
                std::cout << "       alphas_1[t-1][0] " << alphas_1[t-1][0] << " logA[0][k] " << logA[0][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                std::cout << "       alphas_1[t-1][1] " << alphas_1[t-1][1] << " logA[1][k] " << logA[1][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                std::cout << "       alphas_1[t-1][2] " << alphas_1[t-1][2] << " logA[2][k] " << logA[2][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                std::cout << "       alphas_1[t-1][3] " << alphas_1[t-1][3] << " logA[3][k] " << logA[3][k] << " this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                return false;
            }
        }
    }
    return true;
}


// Backward-algorithm: log-space
template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::iBackward(String<String<long double> > &betas_1, unsigned s, unsigned i, String<String<long double> > &logA, AppOptions &options)
{
    unsigned T = this->setObs[s][i].length();
    // for t = T
    for (unsigned k = 0; k < this->K; ++k)
       betas_1[this->setObs[s][i].length() - 1][k] = log(1.0);
    
    // for t = 2:T
    for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
    {
        for (unsigned k = 0; k < this->K; ++k)
        {
            // sum over following states
            long double f1 = betas_1[t+1][0] + logA[k][0] + this->eProbs[s][i][t+1][0];
            long double f2 = betas_1[t+1][1] + logA[k][1] + this->eProbs[s][i][t+1][1];
            long double f3 = betas_1[t+1][2] + logA[k][2] + this->eProbs[s][i][t+1][2];
            long double f4 = betas_1[t+1][3] + logA[k][3] + this->eProbs[s][i][t+1][3];

            betas_1[t][k] = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);

            if (std::isinf(betas_1[t][k]))
            {
                std::cout << "ERROR: betas_1[" << t << "][" << k << "] is " << betas_1[t][k] << std::endl;
                return false;
            }
        }
    }
    return true;
}


// for log-space
// interval-wise to avoid storing alpha_1 and beta_1 values for whole genome
// TODO learn 2-> 2/3 only above threshold !?
template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::computeStatePosteriorsFBupdateTrans(AppOptions &options)
{
    String<String<long double> > logA = this->transMatrix;
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            logA[k_1][k_2] = log(this->transMatrix[k_1][k_2]);
 
    String<String<long double> > p;
    resize(p, this->K, Exact());
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
    {
        SEQAN_OMP_PRAGMA(critical)
        resize(p[k_1], this->K, Exact());
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            p[k_1][k_2] = 0.0;
    }
    long double p_2_2 = 0.0;     // for separate learning of trans. prob from '2' -> '2'
    long double p_2_3 = 0.0;     // for separate learning of trans. prob from '2' -> '3'               

    for (unsigned s = 0; s < 2; ++s)
    {
        bool stop = false;
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1))
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            unsigned T = setObs[s][i].length();
            // forward probabilities
            String<String<long double> > alphas_1;
            resize(alphas_1, T, Exact());
            for (unsigned t = 0; t < T; ++t)
            {
                resize(alphas_1[t], this->K, Exact());
            }
            if (!iForward(alphas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }

            // backward probabilities  
            String<String<long double> > betas_1;
            resize(betas_1, T, Exact());
            for (unsigned t = 0; t < T; ++t)
                resize(betas_1[t], this->K, Exact());

            if (!iBackward(betas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }
           
            // compute state posterior probabilities
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
            {
                long double f1 = alphas_1[t][0] + betas_1[t][0];
                long double f2 = alphas_1[t][1] + betas_1[t][1];
                long double f3 = alphas_1[t][2] + betas_1[t][2];
                long double f4 = alphas_1[t][3] + betas_1[t][3];

                long double norm = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);

                for (unsigned k = 0; k < this->K; ++k)
                {
                    this->statePosteriors[s][k][i][t] = myExp(alphas_1[t][k] + betas_1[t][k] - norm);     // store not in log-space!
 
                    if (std::isnan(this->statePosteriors[s][k][i][t]) || std::isinf(this->statePosteriors[s][k][i][t]) || 
                        this->statePosteriors[s][k][i][t] < 0.0 || this->statePosteriors[s][k][i][t] > 1.0) 
                    {
                        std::cout << "ERROR: state posterior probability is " << this->statePosteriors[s][k][i][t] << "." << std::endl;
                        std::cout << "       s: " << s << " i: " << i << " t: " << t << " k:" << k << std::endl;
                        std::cout << "       alphas_1[t][k]: " << alphas_1[t][k] << " betas_1[t][k]: " << betas_1[t][k] << " norm: " << norm << std::endl;
                        stop = true;
                        continue;
                    }
                }
            }

            // update initial probabilities
            for (unsigned k = 0; k < this->K; ++k)
                this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0]; 

            // compute xi values for interval in preparation for new trans. probs
            String<String<long double> > xis;
            resize(xis, this->K, Exact());
            String<String<long double> > p_i;
            resize(p_i, this->K, Exact());
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1)   
            {
                resize(xis[k_1], this->K, 0.0, Exact());
                resize(p_i[k_1], this->K, 0.0, Exact());
            }
            long double p_2_2_i = 0.0;
            long double p_2_3_i = 0.0;
            // 
            for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
            {
                long double norm = std::numeric_limits<long double>::quiet_NaN();
                for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
                {
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                    {
                        xis[k_1][k_2] = alphas_1[t-1][k_1] + logA[k_1][k_2] + this->eProbs[s][i][t][k_2] + betas_1[t][k_2];
                        norm = get_logSumExp(norm, xis[k_1][k_2], options.lookUp);
                    }
                }
                for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        p_i[k_1][k_2] += myExp(xis[k_1][k_2] - norm);

                // learn p[2->2/3] for region over nThresholdForP
                if (options.nThresholdForTransP > 0 && setObs[s][i].nEstimates[t] >= options.nThresholdForTransP)
                {
                    p_2_2_i += myExp(xis[2][2] - norm);
                    p_2_3_i += myExp(xis[2][3] - norm);
                }
            }
            // add to global sum
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1) 
                for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                    SEQAN_OMP_PRAGMA(critical)
                        p[k_1][k_2] += p_i[k_1][k_2];
            
            SEQAN_OMP_PRAGMA(critical)
                p_2_2 += p_2_2_i;
            SEQAN_OMP_PRAGMA(critical)
                p_2_3 += p_2_3_i;
        }
        if (stop) return false;
    }

    // update transition matrix
    String<String<long double> > A = this->transMatrix;
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
    {
        long double denumerator = 0.0;
        for (unsigned k_3 = 0; k_3 < this->K; ++k_3)
            denumerator += p[k_1][k_3]; 

        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
        {
            A[k_1][k_2] = p[k_1][k_2] / denumerator;
            if (A[k_1][k_2] <= 0.0) A[k_1][k_2] = DBL_MIN;          // make sure not getting zero
        }
    }
    // Fix p[2->2/3] using only trans. probs. for region over nThresholdForP, while keeping sum of p[2->2] and p[2->3] constant 
    if (options.nThresholdForTransP > 0)
    {
        long double sum_2_23 = A[2][2] + A[2][3];
        A[2][2] = sum_2_23 * p_2_2/(p_2_2 + p_2_3);
        A[2][3] = sum_2_23 * p_2_3/(p_2_2 + p_2_3);
    }
    // keep transProb of '2' -> '3' on min. value
    if (A[2][3] < options.minTransProbCS)
    {
        A[2][3] = options.minTransProbCS;

        if (A[3][3] < options.minTransProbCS) A[3][3] = options.minTransProbCS;
        std::cout << "NOTE: Prevented transition probability '2' -> '3' from dropping below min. value of " << options.minTransProbCS << ". Set for transitions '2' -> '3' (and if necessary also for '3'->'3') to " << options.minTransProbCS << "." << std::endl;
    }
    this->transMatrix = A;
    return true;
}


// without updating transition probabilities: log space
template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::computeStatePosteriorsFB(AppOptions &options)
{
    String<String<long double> > logA = this->transMatrix;
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            logA[k_1][k_2] = log(this->transMatrix[k_1][k_2]);

    for (unsigned s = 0; s < 2; ++s)
    {
        bool stop = false;
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            unsigned T = setObs[s][i].length();
            // forward probabilities
            String<String<long double> > alphas_1;
            resize(alphas_1, T, Exact());
            for (unsigned t = 0; t < T; ++t)
            {
                resize(alphas_1[t], this->K, Exact());
            }
            if (!iForward(alphas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }

            // backward probabilities
            String<String<long double> > betas_1;
            resize(betas_1, T, Exact());
            for (unsigned t = 0; t < T; ++t)
                resize(betas_1[t], this->K, Exact());
            if (!iBackward(betas_1, s, i, logA, options))
            {
                stop = true;
                continue;
            }
 
            // compute state posterior probabilities (in log-space)
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
            {
                long double f1 = alphas_1[t][0] + betas_1[t][0];
                long double f2 = alphas_1[t][1] + betas_1[t][1];
                long double f3 = alphas_1[t][2] + betas_1[t][2];
                long double f4 = alphas_1[t][3] + betas_1[t][3];

                long double norm = get_logSumExp_states(f1, f2, f3, f4, options.lookUp);

                for (unsigned k = 0; k < this->K; ++k)
                {
                    this->statePosteriors[s][k][i][t] = myExp(alphas_1[t][k] + betas_1[t][k] - norm);     // store not in log-space!
                    if (std::isnan(this->statePosteriors[s][k][i][t])) std::cout << "ERROR: statePosterior is nan! " << std::endl;
                }
            }

            // update init probs
            for (unsigned k = 0; k < this->K; ++k)
                this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];

        }
        if (stop) return false;
    }
    return true;
}



bool updateDensityParams2(String<String<String<double> > > &statePosteriors1, String<String<String<double> > > &statePosteriors2, String<String<Observations> > &setObs, 
                          GAMMA &gamma1, GAMMA &gamma2, 
                          unsigned &iter, unsigned &trial,
                          AppOptions &options)
{
    if (!gamma1.updateThetaAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
        return false;

    if (!gamma2.updateThetaAndK(statePosteriors2, setObs, options.g2_kMin, options.g2_kMax, options))         // make sure g1k <= g2k
        return false;

    // make sure gamma1.mu < gamma2.mu   
    checkOrderG1G2(gamma1, gamma2, options);
    return true;
}

bool updateDensityParams2(String<String<String<double> > > &statePosteriors1, String<String<String<double> > > &statePosteriors2, String<String<Observations> > &setObs, 
                          GAMMA_REG &gamma1, GAMMA_REG &gamma2, 
                          unsigned &iter, unsigned &trial,
                          AppOptions &options)
{
    if (!gamma1.updateRegCoeffsAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
        return false;

    double g2_kMin = options.g2_kMin;
    if (options.g1_k_le_g2_k)
        g2_kMin = std::max(gamma1.k, options.g2_kMin);

    if (!gamma2.updateRegCoeffsAndK(statePosteriors2, setObs, g2_kMin, options.g2_kMax, options))
        return false;

    // make sure gamma1.mu < gamma2.mu    
    checkOrderG1G2(gamma1, gamma2, iter, trial, options);
    return true;
}


template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::updateDensityParams(TGAMMA &gamma1, TGAMMA &gamma2, unsigned &iter, unsigned &trial, AppOptions &options)   
{
    String<String<String<double> > > statePosteriors1;
    String<String<String<double> > > statePosteriors2;
    resize(statePosteriors1, 2, Exact());
    resize(statePosteriors2, 2, Exact());
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors1[s], length(this->statePosteriors[s][0]), Exact());
        resize(statePosteriors2[s], length(this->statePosteriors[s][0]), Exact());
        for (unsigned i = 0; i < length(this->statePosteriors[s][0]); ++i)
        {
            resize(statePosteriors1[s][i], length(this->statePosteriors[s][0][i]), Exact());
            resize(statePosteriors2[s][i], length(this->statePosteriors[s][0][i]), Exact());
            for (unsigned t = 0; t < length(this->statePosteriors[s][0][i]); ++t)
            {
                statePosteriors1[s][i][t] = this->statePosteriors[s][0][i][t] + this->statePosteriors[s][1][i][t];
                statePosteriors2[s][i][t] = this->statePosteriors[s][2][i][t] + this->statePosteriors[s][3][i][t];
            }
        }
    }

    updateDensityParams2(statePosteriors1, statePosteriors2, this->setObs, gamma1, gamma2, iter, trial, options); 

    return true;
}


template<typename TGAMMA, typename TBIN>
bool HMM<TGAMMA, TBIN>::updateDensityParams(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &options)   
{
    String<String<String<double> > > statePosteriors1;
    String<String<String<double> > > statePosteriors2;
    resize(statePosteriors1, 2, Exact());
    resize(statePosteriors2, 2, Exact());
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors1[s], length(this->statePosteriors[s][0]), Exact());
        resize(statePosteriors2[s], length(this->statePosteriors[s][0]), Exact());
        for (unsigned i = 0; i < length(this->statePosteriors[s][0]); ++i)
        {
            resize(statePosteriors1[s][i], length(this->statePosteriors[s][0][i]), Exact());
            resize(statePosteriors2[s][i], length(this->statePosteriors[s][0][i]), Exact());
            for (unsigned t = 0; t < length(this->statePosteriors[s][0][i]); ++t)
            {
                statePosteriors1[s][i][t] = this->statePosteriors[s][2][i][t];
                statePosteriors2[s][i][t] = this->statePosteriors[s][3][i][t];
            }
        }
    }

    // truncation counts
    modelParams.bin1.updateP(statePosteriors1, this->setObs, options); 
    modelParams.bin2.updateP(statePosteriors2, this->setObs, options);

    // make sure bin1.p < bin2.p   
    checkOrderBin1Bin2(modelParams.bin1, modelParams.bin2);

    return true;
}



// Baum-Welch
// in log-space (using log-sum-exp trick)
template<typename TGAMMA, typename TBIN> 
bool HMM<TGAMMA, TBIN>::baumWelch(ModelParams<TGAMMA, TBIN> &modelParams, CharString learnTag, AppOptions &options)
{
    TGAMMA prev_gamma1 = modelParams.gamma1;
    TGAMMA prev_gamma2 = modelParams.gamma2;
    TBIN prev_bin1 = modelParams.bin1;
    TBIN prev_bin2 = modelParams.bin2;
    unsigned trial = 0;
    for (unsigned iter = 0; iter < options.maxIter_bw; ++iter)
    {
        std::cout << ".. " << iter << "th iteration " << std::endl;
        std::cout << "                        computeEmissionProbs() " << std::endl;
        if (!computeEmissionProbs(modelParams, true, options) )
        {
            std::cerr << "ERROR: Could not compute emission probabilities! " << std::endl;
            return false;
        }
        std::cout << "                        computeStatePosteriorsFB() " << std::endl;
        if (!computeStatePosteriorsFBupdateTrans(options))
        {
            std::cerr << "ERROR: Could not compute forward-backward algorithm! " << std::endl;
            return false;
        }
        
        std::cout << "                        updateDensityParams() " << std::endl;

        if (learnTag == "LEARN_BINOMIAL")
        {
            if (!updateDensityParams(modelParams, options))
            {
                std::cerr << "ERROR: Could not update parameters! " << std::endl;
                return false;
            }
        }
        else
        {
            if (!updateDensityParams(modelParams.gamma1, modelParams.gamma2, iter, trial, options))
            {
                std::cerr << "ERROR: Could not update parameters! " << std::endl;
                return false;
            }
            if (trial > 10)
            {
                std::cerr << "ERROR: Could not learn gamma parameters, exceeded max. number of reseedings! " << std::endl;
                return false;
            }

        }
        
        if (learnTag == "LEARN_GAMMA" && checkConvergence(modelParams.gamma1, prev_gamma1, options) && checkConvergence(modelParams.gamma2, prev_gamma2, options) )             
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        else if (learnTag != "LEARN_GAMMA" && checkConvergence(modelParams.bin1, prev_bin1, options) && checkConvergence(modelParams.bin2, prev_bin2, options) )             
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        prev_gamma1 = modelParams.gamma1;
        prev_gamma2 = modelParams.gamma2;
        prev_bin1 = modelParams.bin1;
        prev_bin2 = modelParams.bin2;

        myPrint(modelParams.gamma1);
        myPrint(modelParams.gamma2);
        
        std::cout << "*** Transition probabilitites ***" << std::endl;
        for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        {
            std::cout << "    " << k_1 << ": ";
            for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                std::cout << this->transMatrix[k_1][k_2] << "  ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
        if (learnTag != "LEARN_GAMMA")
         {
            myPrint(modelParams.bin1);
            myPrint(modelParams.bin2);
        }
    }
    return true;
}


template<typename TGAMMA, typename TBIN> 
bool HMM<TGAMMA, TBIN>::applyParameters(ModelParams<TGAMMA, TBIN> &modelParams, AppOptions &options)
{
    if (!computeEmissionProbs(modelParams, false, options))
    {
        std::cerr << "ERROR: Could not compute emission probabilities! " << std::endl;
        return false;
    }
    if (!computeStatePosteriorsFB(options))
    {
        std::cerr << "ERROR: Could not compute forward-backward algorithm! " << std::endl;
        return false;
    }
    return true;
}



template<typename TGAMMA, typename TBIN>
void HMM<TGAMMA, TBIN>::posteriorDecoding(String<String<String<__uint8> > > &states)
{ 
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(states[s], length(this->setObs[s]), Exact());
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            if (!this->setObs[s][i].discard)
            {
                resize(states[s][i], this->setObs[s][i].length(), Exact());
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                {
                    long double max_p = 0.0;
                    unsigned max_k = 0;
                    for (unsigned k = 0; k < this->K; ++k)
                    {
                        if (this->statePosteriors[s][k][i][t] > max_p)
                        {
                            max_p = this->statePosteriors[s][k][i][t];
                            max_k = k;
                        }
                    }
                    states[s][i][t] = max_k;
                }
            }
        }
    }
}


void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, 
                         String<String<Observations> > &setObs, 
                         String<Data> &data_replicates,
                         String<ModelParams<GAMMA_REG, ZTBIN> > &modelParams)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data_replicates[0].setObs[s]); ++i)
        {
            if (!data_replicates[0].setObs[s][i].discard)
            {
                for (unsigned t = 0; t < data_replicates[0].setObs[s][i].length(); ++t)
                {
                    if (states[s][i][t] >= 2)
                    {
                        double x1 = data_replicates[0].setObs[s][i].rpkms[t];   // same for all replicates
                        bool rm = true;
                        // rm if not > gamma1 expected value in all replicates
                        for (unsigned rep = 0; rep < length(data_replicates); ++rep)
                        {
                            double gamma1_pred = exp(modelParams[rep].gamma1.b0 + modelParams[rep].gamma1.b1 * x1);
                            if (data_replicates[rep].setObs[s][i].kdes[t] > gamma1_pred) rm = false;
                        }
                        if (rm) states[s][i][t] -= 2;
                    }
                }
            }
        }
    }
}

void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, 
                         String<String<Observations> > &setObs, 
                         String<Data> &data_replicates,
                         String<ModelParams<GAMMA_REG, ZTBIN_REG> > &modelParams)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data_replicates[0].setObs[s]); ++i)
        {
            if (!data_replicates[0].setObs[s][i].discard)
            {
                for (unsigned t = 0; t < data_replicates[0].setObs[s][i].length(); ++t)
                {
                    if (states[s][i][t] >= 2)
                    {
                        double x1 = data_replicates[0].setObs[s][i].rpkms[t];   // same for all replicates
                        bool rm = true;
                        // rm if not > gamma1 expected value in all replicates
                        for (unsigned rep = 0; rep < length(data_replicates); ++rep)
                        {
                            double gamma1_pred = exp(modelParams[rep].gamma1.b0 + modelParams[rep].gamma1.b1 * x1);
                            if (data_replicates[rep].setObs[s][i].kdes[t] > gamma1_pred) rm = false;
                        }
                        if (rm) states[s][i][t] -= 2;
                    }
                }
            }
        }
    }
}

// do nothing
void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, 
                         String<String<Observations> > &setObs, 
                         String<Data> &data_replicates, 
                         String<ModelParams<GAMMA, ZTBIN> > &modelParams){}

void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, 
                         String<String<Observations> > &setObs, 
                         String<Data> &data_replicates, 
                         String<ModelParams<GAMMA, ZTBIN_REG> > &modelParams){}


// for GLM with input signal: 
// when using free gamma shapes, i.e. gamma1.k can be > gamma2.k
// make sure sites with fragment coverage (KDE) below gamma1.mean are classified as 'non-enriched'
template<typename TGAMMA, typename TBIN>
void HMM<TGAMMA, TBIN>::rmBoarderArtifacts(String<String<String<__uint8> > > &states, 
                                           String<Data> &data_replicates, 
                                           String<ModelParams<TGAMMA, TBIN> > &modelParams)
{
    rmBoarderArtifacts2(states, this->setObs, data_replicates, modelParams);
}


void writeStates(String<BedRecord<Bed6> > &bedRecords_sites,
                 Data &data,
                 FragmentStore<> &store, 
                 unsigned contigId,
                 AppOptions &options)          
{ 
    long double min_val;
    if (!options.useHighPrecision)
        min_val = DBL_MIN;
    else
        min_val = LDBL_MIN;

    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)    // data.states[s]
        {
            // note: could skip discarded intervals here ...
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)        // length(data.states[s][i])
            {
                if (options.outputAll && data.setObs[s][i].truncCounts[t] >= 1 && !data.setObs[s][i].discard)
                {
                    BedRecord<Bed6> record;

                    record.ref = store.contigNameStore[contigId];
                    if (s == 0)         // '+'-strand; crosslink sites (not truncation site)
                    {
                        if (!options.crosslinkAtTruncSite)  // default
                            record.beginPos = t + data.setPos[s][i] - 1;
                        else
                            record.beginPos = t + data.setPos[s][i];

                        record.endPos = record.beginPos + 1;
                    }
                    else                 // '-'-strand;
                    {
                        if (!options.crosslinkAtTruncSite) 
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                        else
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1;

                        record.endPos = record.beginPos + 1;
                    }

                    std::stringstream ss;
                    ss << (int)data.states[s][i][t];
                    record.name = ss.str();
                    ss.str("");  
                    ss.clear();  

                    if (options.score_type == 0)
                    {
                        // log posterior prob. ratio score
                        long double secondBest = 0.0;
                        for (unsigned k = 0; k < 4; ++k)
                        {
                            if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                                secondBest = data.statePosteriors[s][k][i][t];
                        }      
                        ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max((long double)secondBest, min_val));
                    }
                    else if (options.score_type == 1)
                    {
                        // log(3/2)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][2][i][t], min_val));
                    }
                    else if (options.score_type == 2)
                    {
                        // log(3/1)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][1][i][t], min_val));
                    }
                    else if (options.score_type == 3)
                    {
                        // log(3/0)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][0][i][t], min_val));
                    }                    
                    else if (options.score_type == 4)
                    {
                        // log(3/(1+2))  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][1][i][t] + (long double)data.statePosteriors[s][2][i][t], min_val));
                    }
                    else if (options.score_type == 5)
                    {
                        // 3 
                        ss << (double)data.statePosteriors[s][3][i][t];
                    }
                    else if (options.score_type == 6)
                    {
                        // log(enriched/non-enriched) + log(crosslinked/non-crosslinked)  
                        ss << ((double)log((data.statePosteriors[s][2][i][t] + data.statePosteriors[s][3][i][t])/std::max((long double)data.statePosteriors[s][0][i][t] + (long double)data.statePosteriors[s][1][i][t], min_val)) + (double)log((data.statePosteriors[s][1][i][t] + data.statePosteriors[s][3][i][t])/std::max((long double)data.statePosteriors[s][0][i][t] + (long double)data.statePosteriors[s][2][i][t], min_val))) ;
                    }

                    record.score = ss.str();
                    ss.str("");  
                    ss.clear();  
                    if (s == 0)
                        record.strand = '+';
                    else
                        record.strand = '-';

                    ss << 0;
                    ss << ";";
                    ss << (int)data.setObs[s][i].truncCounts[t];
                    ss << ";";
                    ss << (int)data.setObs[s][i].nEstimates[t];
                    ss << ";";
                    ss << (double)data.setObs[s][i].kdes[t];
                    ss << ";";

                    ss << (double)data.statePosteriors[s][3][i][t];
                    ss << ";"; 
                    if (options.useCov_RPKM)
                        ss << (double)data.setObs[s][i].rpkms[t];
                    else
                        ss << 0.0;
                    ss << ";";
                    ss << (double)log((data.statePosteriors[s][2][i][t] + data.statePosteriors[s][3][i][t])/(data.statePosteriors[s][0][i][t] + data.statePosteriors[s][1][i][t]));
                    ss << ";";

                    record.data = ss.str();
                    ss.str("");  
                    ss.clear();  

                    appendValue(bedRecords_sites, record);
                }
                else if (data.setObs[s][i].discard && options.outputAll && data.setObs[s][i].truncCounts[t] >= 1)  // discarded interval
                {
                    BedRecord<Bed6> record;

                    record.ref = store.contigNameStore[contigId];
                    if (s == 0)         // '+'-strand; crosslink sites (not truncation site)
                    {
                        if (!options.crosslinkAtTruncSite)  // default
                            record.beginPos = t + data.setPos[s][i] - 1;
                        else
                            record.beginPos = t + data.setPos[s][i];

                        record.endPos = record.beginPos + 1;
                    }
                    else                 // '-'-strand;
                    {
                        if (!options.crosslinkAtTruncSite) 
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                        else
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1;

                        record.endPos = record.beginPos + 1;
                    }

                    std::stringstream ss;
                    ss << (int)0;           // assign 'non-enriched + non-crosslink' 
                    record.name = ss.str();
                    ss.str("");  
                    ss.clear();  

                    // log posterior prob. ratio score                   
                    ss << "NA";

                    record.score = ss.str();
                    ss.str("");  
                    ss.clear();  
                    if (s == 0)
                        record.strand = '+';
                    else
                        record.strand = '-';

                    ss << 0;
                    ss << ";";
                    ss << (int)data.setObs[s][i].truncCounts[t];
                    ss << ";";
                    ss << (int)data.setObs[s][i].nEstimates[t];
                    ss << ";";
                    ss << (double)data.setObs[s][i].kdes[t];
                    ss << ";";

                    ss << "NA";
                    ss << ";"; 
                    if (options.useCov_RPKM)
                        ss << (double)data.setObs[s][i].rpkms[t];
                    else
                        ss << 0.0;
                    ss << ";";
                    ss << "NA";
                    ss << ";";

                    record.data = ss.str();
                    ss.str("");  
                    ss.clear();  

                    appendValue(bedRecords_sites, record);
                }
                else if (!data.setObs[s][i].discard && data.states[s][i][t] == 3)
                {
                    BedRecord<Bed6> record;

                    record.ref = store.contigNameStore[contigId];
                    if (s == 0)         // crosslink sites (not truncation site)
                    {
                        if (!options.crosslinkAtTruncSite)
                            record.beginPos = t + data.setPos[s][i] - 1;
                        else
                            record.beginPos = t + data.setPos[s][i];

                        record.endPos = record.beginPos + 1;
                    }
                    else
                    {
                        if (!options.crosslinkAtTruncSite)
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                        else
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1;

                        record.endPos = record.beginPos + 1;
                    }

                    std::stringstream ss;
                    ss << (int)data.states[s][i][t];
                    record.name = ss.str();
                    ss.str("");  
                    ss.clear();  

                    if (options.score_type == 0)
                    {
                        // log posterior prob. ratio score
                        long double secondBest = 0.0;
                        for (unsigned k = 0; k < 4; ++k)
                        {
                            if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                                secondBest = data.statePosteriors[s][k][i][t];
                        }      
                        ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) );
                    }
                    else if (options.score_type == 1)
                    {
                        // log(3/2)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][2][i][t], min_val));
                    }
                    else if (options.score_type == 2)
                    {
                        // log(3/1)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][1][i][t], min_val));
                    }
                    else if (options.score_type == 3)
                    {
                        // log(3/0)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][0][i][t], min_val));
                    }                    
                    else if (options.score_type == 4)
                    {
                        // log(3/(1+2))  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/std::max((long double)data.statePosteriors[s][1][i][t] + (long double)data.statePosteriors[s][2][i][t], min_val));
                    }
                    else if (options.score_type == 5)
                    {
                        // 3 
                        ss << (double)data.statePosteriors[s][3][i][t];
                    }
                    else if (options.score_type == 6)
                    {
                        // log(enriched/non-enriched) + log(crosslinked/non-crosslinked)  
                        ss << ((double)log((data.statePosteriors[s][2][i][t] + data.statePosteriors[s][3][i][t])/std::max((long double)data.statePosteriors[s][0][i][t] + (long double)data.statePosteriors[s][1][i][t], min_val)) + (double)log((data.statePosteriors[s][1][i][t] + data.statePosteriors[s][3][i][t])/std::max((long double)data.statePosteriors[s][0][i][t] + (long double)data.statePosteriors[s][2][i][t], min_val))) ;
                    }

                    record.score = ss.str();
                    ss.str("");  
                    ss.clear();  
                    if (s == 0)
                        record.strand = '+';
                    else
                        record.strand = '-';

                    appendValue(bedRecords_sites, record);
                }               
            }
        }
    }
}


void writeRegions(String<BedRecord<Bed6> > &bedRecords_regions,
                 Data &data,
                 FragmentStore<> &store, 
                 unsigned contigId,
                 AppOptions &options)          
{ 
    long double min_val;
    if (!options.useHighPrecision)
        min_val = DBL_MIN;
    else
        min_val = LDBL_MIN;

    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.states[s]); ++i)
        {
            for (unsigned t = 0; t < length(data.states[s][i]); ++t)
            {
                if (!data.setObs[s][i].discard && data.states[s][i][t] == 3)
                {
                    BedRecord<Bed6> record;
                    record.ref = store.contigNameStore[contigId];
                    if (s == 0)         // crosslink sites (not truncation site)
                    {
		                if (!options.crosslinkAtTruncSite)
			                record.beginPos = t + data.setPos[s][i] - 1;
                        else
			                record.beginPos = t + data.setPos[s][i];
                        
                        record.endPos = record.beginPos + 1;
                    }
                    else
                    {
		                if (!options.crosslinkAtTruncSite)
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                        else
                            record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1;

                        record.endPos = record.beginPos + 1;
                    }

                    std::stringstream ss;

                    // TODO use different scores
                    // log posterior prob. ratio score
                    long double secondBest = 0.0;
                    for (unsigned k = 0; k < 4; ++k)
                    {
                        if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                            secondBest = data.statePosteriors[s][k][i][t];
                    }
                    ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) );

                    record.score = ss.str();
                    ss.str("");  
                    ss.clear();  
                    if (s == 0)
                        record.strand = '+';
                    else
                        record.strand = '-';

                    unsigned prev_cs = t;
                    long double scoresSum = (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) );
                    std::stringstream ss_indivScores;
                    ss_indivScores << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) ) << ';';
                    while ((t+1) < length(data.states[s][i]) && (t+1-prev_cs) <= options.distMerge)
                    {
                        ++t;
                        if (data.states[s][i][t] == 3)
                        {
                            if (s == 0)         // crosslink sites (not truncation site)
                            {
                                if (!options.crosslinkAtTruncSite)
                                    record.endPos = t + data.setPos[s][i];  
                                else
                                    record.endPos = t + data.setPos[s][i] + 1;  
                            }
                            else
                            {
                                if (!options.crosslinkAtTruncSite)                                
                                    record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]);
                                else
                                    record.beginPos = length(store.contigStore[contigId].seq) - (t + data.setPos[s][i]) - 1;
                            }

                            // log posterior prob. ratio score
                            long double secondBest = 0.0;
                            for (unsigned k = 0; k < 4; ++k)
                            {
                                if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                                    secondBest = data.statePosteriors[s][k][i][t];
                            }                   

                            scoresSum += (long double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) );
                            ss_indivScores << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) ) << ';';
                            prev_cs = t;
                        }
                    }
                    ss << scoresSum;
                    record.score = ss.str();
                    ss.str("");  
                    ss.clear(); 
                    record.name = ss_indivScores.str();
                    ss_indivScores.str("");  
                    ss_indivScores.clear(); 
                    
                    appendValue(bedRecords_regions, record);
                }      
            }
        }
    }
}

  
template<typename TGAMMA, typename TBIN>
void myPrint(HMM<TGAMMA, TBIN> &hmm)
{
    std::cout << "*** Transition probabilitites ***" << std::endl;
    for (unsigned k_1 = 0; k_1 < hmm.K; ++k_1)
    {
        std::cout << "    " << k_1 << ": ";
        for (unsigned k_2 = 0; k_2 < hmm.K; ++k_2)
            std::cout << hmm.transMatrix[k_1][k_2] << "  ";
        std::cout << std::endl;
    }
}


template<typename TOut>
void printParams(TOut &out, String<String<long double> > &transMatrix)
{
    out << "Transition probabilities:" << std::endl;
    for (unsigned k_1 = 0; k_1 < 4; ++k_1)
    {
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
            out << transMatrix[k_1][k_2] << "\t";
        out << std::endl;
    }
}

#endif
