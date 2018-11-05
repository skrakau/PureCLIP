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


template <typename TGAMMA, typename TBIN, typename TDOUBLE>
class HMM {     

public:

    __uint8                                 K;                  // no. of sates
    String<String<String<long double> > >   initProbs;          // intital probabilities

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
//     HMM<TGAMMA, TBIN, TDOUBLE>(const HMM<TGAMMA, TBIN, TDOUBLE> &hmm2) {x = p2.x; y = p2.y; } 

    HMM<TGAMMA, TBIN, TDOUBLE>();
    ~HMM<TGAMMA, TBIN, TDOUBLE>();
    void setInitProbs(String<double> &probs);
    bool computeEmissionProbs(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, bool learning, AppOptions &options);
    bool iForward(String<String<TDOUBLE> > &alphas_1, String<String<TDOUBLE> > &alphas_2, unsigned s, unsigned i, AppOptions &options);
    bool iBackward(String<String<TDOUBLE> > &betas_2, String<String<TDOUBLE> > &alphas_1, unsigned s, unsigned i);
    bool computeStatePosteriorsFB(AppOptions &options);
    bool computeStatePosteriorsFBupdateTrans(AppOptions &options);
    bool updateTransProbsAndPostDecode(AppOptions &options);
    bool updateDensityParams(TGAMMA &d1, TGAMMA &d2, AppOptions &options);
    bool updateDensityParams(TGAMMA /*&d1*/, TGAMMA /*&d2*/, TBIN &bin1, TBIN &bin2, AppOptions &options);
    bool baumWelch(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, CharString learnTag, AppOptions &options);
    bool applyParameters(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, AppOptions &/*options*/);
    long double viterbi(String<String<String<__uint8> > > &states);
    long double viterbi_log(String<String<String<__uint8> > > &states);
    void posteriorDecoding(String<String<String<__uint8> > > &states);
    void rmBoarderArtifacts(String<String<String<__uint8> > > &states, TGAMMA &g1);

    // for each F/R,interval,t, state ....
    String<String<String<String<TDOUBLE> > > > eProbs;           // emission/observation probabilities  P(Y_t | S_t) -> precompute for each t given Y_t = (C_t, T_t) !!!
    String<String<String<String<TDOUBLE> > > > statePosteriors;  // for each k: for each covered interval string of posteriors
};


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
HMM<TGAMMA, TBIN, TDOUBLE>::~HMM<TGAMMA, TBIN, TDOUBLE>()
{
    clear(this->eProbs);
    clear(this->statePosteriors);
    clear(this->initProbs);
    clear(this->transMatrix);
   // do not touch observations
}

 
template<typename TEProbs, typename TSetObs, typename TDOUBLE>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA2<TDOUBLE> &d1, GAMMA2<TDOUBLE> &d2, ZTBIN<TDOUBLE> &bin1, ZTBIN<TDOUBLE> &bin2, unsigned t, AppOptions &options)
{
    long double g1_d = 1.0;
    long double g2_d = 0.0;
    if (setObs.kdes[t] >= d1.tp) 
    {
        g1_d = d1.getDensity(setObs.kdes[t]);
        g2_d = d2.getDensity(setObs.kdes[t]); 
    }
    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);

        if (options.keep_unlikelyNkratios)
        {
            // for high k/N ratios, if eProbs -> 0.0, do not discard interval! 
            // Force states is at boarders of binomial distributions using pseudo-eProbs! 
            double nk_ratio = (setObs.nEstimates[t] > setObs.truncCounts[t]) ? ((double)setObs.truncCounts[t]/(double)setObs.nEstimates[t]) : ((double)setObs.truncCounts[t]/setObs.truncCounts[t]);
            if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin1_d == 0.0) && 
                    (nk_ratio > bin2.p))
            {
                // ok, if crosslink score -> inf; default score -> log(P('3'|Y)/P('1'|Y)) or  log(1.0/DBL_MIN) 
                bin2_d = 1.0;
            }
            else if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin2_d == 0.0) && 
                    (nk_ratio < bin1.p))
            {
                bin1_d = 1.0;  
            }
        }
        // else if 0s, will be discarded downstream ...
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    // 
    if ((eProbs[0] == 0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0) ||
            (g1_d + g2_d < options.min_eProbSum) || (bin1_d + bin2_d < options.min_eProbSum) ||
            (eProbs[0] + eProbs[1] + eProbs[2] + eProbs[3] < options.min_eProbSum) ||
            (std::isnan(eProbs[0]) || std::isnan(eProbs[1]) || std::isnan(eProbs[2]) || std::isnan(eProbs[3])) ||
            (std::isinf(eProbs[0]) || std::isinf(eProbs[1]) || std::isinf(eProbs[2]) || std::isinf(eProbs[3])))
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
                std::cout << "       emission probability 'non-enriched' gamma: " << g1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << g2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 1.0;
        eProbs[1] = 0.0;
        eProbs[2] = 0.0;
        eProbs[3] = 0.0;
        
        return false;
    }
    return true;
}

template<typename TEProbs, typename TSetObs, typename TDOUBLE>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA2_REG<TDOUBLE> &d1, GAMMA2_REG<TDOUBLE> &d2, ZTBIN<TDOUBLE> &bin1, ZTBIN<TDOUBLE> &bin2, unsigned t, AppOptions &options)
{
    long double x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
    long double d1_pred = exp(d1.b0 + d1.b1 * x);
    long double d2_pred = exp(d2.b0 + d2.b1 * x);

    long double g1_d = 1.0;
    long double g2_d = 0.0;
    if (setObs.kdes[t] >= d1.tp) 
    {
        g1_d = d1.getDensity(setObs.kdes[t], d1_pred, options);
        g2_d = d2.getDensity(setObs.kdes[t], d2_pred, options); 
    }

    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], options);

        if (options.keep_unlikelyNkratios)
        {
            // for high k/N ratios, if eProbs -> 0.0, do not discard interval! 
            // Force states is at boarders of binomial distributions using pseudo-eProbs! 
            double nk_ratio = (setObs.nEstimates[t] > setObs.truncCounts[t]) ? ((double)setObs.truncCounts[t]/(double)setObs.nEstimates[t]) : ((double)setObs.truncCounts[t]/setObs.truncCounts[t]);
            if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin1_d == 0.0) && 
                    (nk_ratio > bin2.p))
            {
                // ok, if crosslink score -> inf; default score -> log(P('3'|Y)/P('1'|Y)) or  log(1.0/DBL_MIN) 
                bin2_d = 1.0;
            }
            else if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin2_d == 0.0) && 
                    (nk_ratio < bin1.p))
            {
                bin1_d = 1.0;  
            }
        }
        // else if 0s, will be discarded downstream ...
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    // 
    if ((eProbs[0] == 0.0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0) || 
            (g1_d + g2_d < options.min_eProbSum) || (bin1_d + bin2_d < options.min_eProbSum) ||
            (eProbs[0] + eProbs[1] + eProbs[2] + eProbs[3] < options.min_eProbSum) ||
            (std::isnan(eProbs[0]) || std::isnan(eProbs[1]) || std::isnan(eProbs[2]) || std::isnan(eProbs[3])) ||
            (std::isinf(eProbs[0]) || std::isinf(eProbs[1]) || std::isinf(eProbs[2]) || std::isinf(eProbs[3])))
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
                std::cout << "       covariate b: " << x << " predicted mean 'non-enriched': " << d1_pred << " predicted mean 'enriched': " << d2_pred << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << g1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << g2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }    
        eProbs[0] = 1.0;
        eProbs[1] = 0.0;
        eProbs[2] = 0.0;
        eProbs[3] = 0.0;

        return false;
    }
    return true;
}

template<typename TEProbs, typename TSetObs, typename TDOUBLE>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA2<TDOUBLE> &d1, GAMMA2<TDOUBLE> &d2, ZTBIN_REG<TDOUBLE> &bin1, ZTBIN_REG<TDOUBLE> &bin2, unsigned t, AppOptions &options)
{
    long double g1_d = 1.0;
    long double g2_d = 0.0;
    if (setObs.kdes[t] >= d1.tp) 
    {
        g1_d = d1.getDensity(setObs.kdes[t]);
        g2_d = d2.getDensity(setObs.kdes[t]); 
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

        if (options.keep_unlikelyNkratios)
        {        
            // for high k/N ratios, if eProbs -> 0.0, do not discard interval! 
            // Force states is at boarders of binomial distributions using pseudo-eProbs! 
            double nk_ratio = (setObs.nEstimates[t] > setObs.truncCounts[t]) ? ((double)setObs.truncCounts[t]/(double)setObs.nEstimates[t]) : ((double)setObs.truncCounts[t]/setObs.truncCounts[t]);
            if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin1_d == 0.0) && 
                    (nk_ratio > bin2_pred))
            {
                // ok, if crosslink score -> inf; default score -> log(P('3'|Y)/P('1'|Y)) or  log(1.0/DBL_MIN) 
                bin2_d = 1.0;
            }
            else if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin2_d == 0.0) && 
                    (nk_ratio < bin1_pred))
            {
                bin1_d = 1.0;  
            }
        }
        // else if 0s, will be discarded downstream ...
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    //
    if ((eProbs[0] == 0.0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0) ||
            (g1_d + g2_d < options.min_eProbSum) || (bin1_d + bin2_d < options.min_eProbSum) || 
            (eProbs[0] + eProbs[1] + eProbs[2] + eProbs[3] < options.min_eProbSum) ||
            (std::isnan(eProbs[0]) || std::isnan(eProbs[1]) || std::isnan(eProbs[2]) || std::isnan(eProbs[3])) ||
            (std::isinf(eProbs[0]) || std::isinf(eProbs[1]) || std::isinf(eProbs[2]) || std::isinf(eProbs[3])) )
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
                std::cout << "       emission probability 'non-enriched' gamma: " << g1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << g2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 1.0;
        eProbs[1] = 0.0;
        eProbs[2] = 0.0;
        eProbs[3] = 0.0;

        return false;
    }
    return true;
}

template<typename TEProbs, typename TSetObs, typename TDOUBLE>
bool computeEProb(TEProbs &eProbs, TSetObs &setObs, GAMMA2_REG<TDOUBLE> &d1, GAMMA2_REG<TDOUBLE> &d2, ZTBIN_REG<TDOUBLE> &bin1, ZTBIN_REG<TDOUBLE> &bin2, unsigned t, AppOptions &options)
{
    long double x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
    long double d1_pred = exp(d1.b0 + d1.b1 * x);
    long double d2_pred = exp(d2.b0 + d2.b1 * x);

    long double g1_d = 1.0;
    long double g2_d = 0.0;
    if (setObs.kdes[t] >= d1.tp) 
    {
        g1_d = d1.getDensity(setObs.kdes[t], d1_pred, options);
        g2_d = d2.getDensity(setObs.kdes[t], d2_pred, options); 
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

        if (options.keep_unlikelyNkratios)
        {
            // for high k/N ratios, if eProbs -> 0.0, do not discard interval! 
            // Force states is at boarders of binomial distributions using pseudo-eProbs! 
            double nk_ratio = (setObs.nEstimates[t] > setObs.truncCounts[t]) ? ((double)setObs.truncCounts[t]/(double)setObs.nEstimates[t]) : ((double)setObs.truncCounts[t]/setObs.truncCounts[t]);
            if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin1_d == 0.0) && 
                    (nk_ratio > bin2_pred))
            {
                // ok, if crosslink score -> inf; default score -> log(P('3'|Y)/P('1'|Y)) or  log(1.0/DBL_MIN) 
                bin2_d = 1.0;
            }
            else if ((g1_d*bin1_d + g1_d*bin2_d + g2_d*bin1_d + g2_d*bin2_d < options.min_eProbSum) && 
                    (bin2_d == 0.0) && 
                    (nk_ratio < bin1_pred))
            {
                bin1_d = 1.0;  
            }
        }
        // else if 0s, will be discarded downstream ...
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    // 
    if ((eProbs[0] == 0.0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0) || 
            (g1_d + g2_d < options.min_eProbSum) || (bin1_d + bin2_d < options.min_eProbSum) ||
            (eProbs[0] + eProbs[1] + eProbs[2] + eProbs[3] < options.min_eProbSum) ||
            (std::isnan(eProbs[0]) || std::isnan(eProbs[1]) || std::isnan(eProbs[2]) || std::isnan(eProbs[3])) ||
            (std::isinf(eProbs[0]) || std::isinf(eProbs[1]) || std::isinf(eProbs[2]) || std::isinf(eProbs[3])))
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
                std::cout << "       covariate b: " << x << " predicted mean 'non-enriched': " << d1_pred << " predicted mean 'enriched': " << d2_pred << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate x: " << setObs.fimoScores[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << g1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << g2_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1_d << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2_d << std::endl;
        }
        eProbs[0] = 1.0;
        eProbs[1] = 0.0;
        eProbs[2] = 0.0;
        eProbs[3] = 0.0;

        return false;
    }
    return true;
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::computeEmissionProbs(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, bool learning, AppOptions &options)
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

                if (!computeEProb(this->eProbs[s][i][t], this->setObs[s][i], d1, d2, bin1, bin2, t, options))
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
                if (!options.useHighPrecision)
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
                if (!options.useHighPrecision)
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



// for one interval only
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::iForward(String<String<TDOUBLE> > &alphas_1, String<String<TDOUBLE> > &alphas_2, unsigned s, unsigned i, AppOptions &options)
{
    // for t = 1
    long double norm = 0.0;
    for (unsigned k = 0; k < this->K; ++k)
    {
        alphas_1[0][k] = this->initProbs[s][i][k] * this->eProbs[s][i][0][k];
        norm += alphas_1[0][k];
    }
    if (norm == 0.0) 
    {
        std::cout << "ERROR: while computing forward values. All alpha values are 0 at t: "<< 0 << "  i: " << i << "." << std::endl;
        if (options.verbosity >= 2)
        {
            for (unsigned k = 0; k < this->K; ++k)
            {
                std::cout << "k: " << k << std::endl;
                std::cout << "this->eProbs[s][i][0][k] " << this->eProbs[s][i][0][k] << std::endl;
            }
        }
        return false;
    }

    for (unsigned k = 0; k < this->K; ++k)
       alphas_2[0][k] = alphas_1[0][k] / norm;

    // for t = 2:T
    for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
    {
        norm = 0.0;
        for (unsigned k = 0; k < this->K; ++k)
        {
            // sum over previous states
            long double sum = 0.0;
            for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                sum += alphas_2[t-1][k_2] * this->transMatrix[k_2][k];
            
            if (sum == 0.0 || std::isnan(sum)) 
            {
                std::cout << "ERROR: while computing forward values. Sum over previous states is " << sum << " at t: "<< t << "  i: " << i << "." << std::endl;
                if (options.verbosity >= 2)
                {
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        std::cout << " k_2: " << k_2 <<  " alphas_2[t-1][k_2]: " << alphas_2[t-1][k_2] << " transMatrix[k_2][k]: " << this->transMatrix[k_2][k] << std::endl; 
                }
                return false;
            }

            // alpha_1
            alphas_1[t][k] = sum * this->eProbs[s][i][t][k];        // - nan?
            norm += alphas_1[t][k];
        }
        
        if (norm == 0.0 || std::isnan(norm)) 
        {
            std::cout << "ERROR: while computing forward values. Sum of alpha values is " << norm << " at t: "<< t << "  i: " << i << "." << std::endl;            
            if (options.verbosity >= 2)
            {
                for (unsigned k = 0; k < this->K; ++k)
                {
                    std::cout << "k: " << k << std::endl;
                    std::cout << "this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                }
            }
            return false;
        }
        // normalize
        for (unsigned k = 0; k < this->K; ++k)
            alphas_2[t][k] = alphas_1[t][k] / norm;   // TODO store scaling coefficients too !?
    }
    return true;
}


// Backward-algorithm
// need alphas_1 for scaling here,
// only betas_2 is needed to compute posterior probs.
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::iBackward(String<String<TDOUBLE> > &betas_2, String<String<TDOUBLE> > &alphas_1, unsigned s, unsigned i)
{
    unsigned T = this->setObs[s][i].length();
    String<String<TDOUBLE> > betas_1;
    resize(betas_1, T, Exact());
    for (unsigned t = 0; t < T; ++t)
        resize(betas_1[t], this->K, Exact());

    // for t = T
    for (unsigned k = 0; k < this->K; ++k)
       betas_1[this->setObs[s][i].length() - 1][k] = 1.0;
    
    long double norm = 0.0;      // use scaling coefficients from alphas here !
    for (unsigned k = 0; k < this->K; ++k)
       norm += alphas_1[this->setObs[s][i].length() - 1][k];

    if (norm == 0.0 || std::isnan(norm)) 
    {
        std::cout << "ERROR: while computing backward values. Sum over alpha values is " << norm << " at t: "<< (this->setObs[s][i].length() - 1) << "  i: " << i << "." << std::endl;
        for (unsigned k = 0; k < this->K; ++k)
            std::cout << "alphas_1[this->setObs[s][i].length() - 1][" << k << "]: " <<  alphas_1[this->setObs[s][i].length() - 1][k] << std::endl;
        return false;
    }

    for (unsigned k = 0; k < this->K; ++k)
    {
       betas_2[this->setObs[s][i].length() - 1][k] = betas_1[this->setObs[s][i].length() - 1][k] / norm;

       if (std::isnan(betas_2[this->setObs[s][i].length() - 1][k]) || std::isinf(betas_2[this->setObs[s][i].length() - 1][k])) 
       {
           std::cout << "ERROR: while computing backward values! betas_2[this->setObs[s][i].length() - 1][k] is " << betas_2[this->setObs[s][i].length() - 1][k] << " at t: "<< (this->setObs[s][i].length() - 1) << "  i: " << i << "." << std::endl;
           std::cout << "betas_1[this->setObs[s][i].length() - 1][k]: " << betas_1[this->setObs[s][i].length() - 1][k] << " norm: " <<  norm << std::endl;
           return false;
       }
    }

    // for t = 2:T
    for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
    {
        norm = 0.0;
        for (unsigned k = 0; k < this->K; ++k)      // precompute ???
            norm += alphas_1[t][k];

        if (norm == 0.0 || std::isnan(norm)) 
        {
            std::cout << "ERROR: while computing backward values. Sum over alpha values is " << norm << " at t: "<< t << "  i: " << i << "." << std::endl;
            for (unsigned k = 0; k < this->K; ++k)
                std::cout << "alphas_1[t][" << k << "]: " <<  alphas_1[t][k] << std::endl;
            return false;
        }

        for (unsigned k = 0; k < this->K; ++k)
        {
            // sum over previous states
            long double sum = 0.0;
            for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                sum += betas_2[t+1][k_2] * this->transMatrix[k][k_2] * this->eProbs[s][i][t+1][k_2];
           
            if (std::isnan(sum))
            {
                std::cout << "ERROR: while computing backward values! sum is nan at t: "<< t << "  i: " << i << " and k: " << k << "." << std::endl;
                for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                    std::cout << "k_2: " << k_2 << " betas_2[t+1][k_2]: " << betas_2[t+1][k_2] << " this->transMatrix[k][k_2]: " <<  this->transMatrix[k][k_2] << " this->eProbs[s][i][t+1][k_2]: " << this->eProbs[s][i][t+1][k_2] << std::endl;
                return false;
            }

            // beta_1
            betas_1[t][k] = sum;
            // beta_2
            betas_2[t][k] = betas_1[t][k] / norm;

            if (std::isnan(betas_2[t][k]) || std::isinf(betas_2[t][k])) 
            {
                std::cout << "ERROR: while computing backward values! betas_2[t][k] is " << betas_2[t][k] << " at t: "<< t << "  i: " << i << "." << std::endl;
                std::cout << "betas_1[t][k]: " << betas_1[t][k] << " norm: " <<  norm << std::endl;
                return false;
            }

        }
    }
    return true;
}


// both for scaling and no-scaling method
/*template<typename TD1, typename TD2, typename TB1, typename TB2>
void HMM<TGAMMA, TBIN, TDOUBLE>::computeStatePosteriors()
{
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
            {
                double sum = 0.0;
                for (unsigned k = 0; k < this->K; ++k)
                    sum += this->alphas_2[s][i][t][k] * this->betas_2[s][i][t][k];

                if (sum == 0.0) 
                {
                    std::cerr << "ERROR: sum == 0 at i: " << i << " t: "<< t << std::endl;
                    for (unsigned k = 0; k < this->K; ++k)
                    {
                        std::cout << "k: " << k << std::endl;
                        std::cout << "this->alphas_2[s][i][t][k]: " << this->alphas_2[s][i][t][k] << " this->betas_2[s][i][t][k]: " << this->betas_2[s][i][t][k] << std::endl;
                    }
                }
                for (unsigned k = 0; k < this->K; ++k)
                    this->statePosteriors[s][k][i][t] = this->alphas_2[s][i][t][k] * this->betas_2[s][i][t][k] / sum;
            }
        }
    }
}*/




// for scaling method
// interval-wise to avoid storing alpha_1, alpha_2 and beta_1, beta_2 values for whole genome
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::computeStatePosteriorsFBupdateTrans(AppOptions &options)
{
    String<String<long double> > A = this->transMatrix;
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
            String<String<TDOUBLE> > alphas_1;
            String<String<TDOUBLE> > alphas_2;
            resize(alphas_1, T, Exact());
            resize(alphas_2, T, Exact());
            for (unsigned t = 0; t < T; ++t)
            {
                resize(alphas_1[t], this->K, Exact());
                resize(alphas_2[t], this->K, Exact());
            } 
            if (!iForward(alphas_1, alphas_2, s, i, options))
            {
                stop = true;
                continue;
            }

            // backward probabilities  
            String<String<TDOUBLE> > betas_2;
            resize(betas_2, T, Exact());
            for (unsigned t = 0; t < T; ++t)
                resize(betas_2[t], this->K, Exact());
            if (!iBackward(betas_2, alphas_1, s, i))
            {
                stop = true;
                continue;
            }
           
            // compute state posterior probabilities
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
            {
                long double sum = 0.0;
                for (unsigned k = 0; k < this->K; ++k)
                    sum += alphas_2[t][k] * betas_2[t][k];

                if (sum == 0.0) 
                {
                    std::cout << "ERROR: while computing state posterior probabilities! Sum of alpha and beta values is 0 at i: " << i << " t: "<< t << "." << std::endl;
                    if (options.verbosity >= 2)
                    {
                        for (unsigned k = 0; k < this->K; ++k)
                        {
                            std::cout << "k: " << k << std::endl;
                            std::cout << "alphas_2[k]: " << alphas_2[t][k] << " betas_2[t][k]: " << betas_2[t][k] << std::endl;
                        }
                    }
                    stop = true;
                    continue;
                }
                for (unsigned k = 0; k < this->K; ++k)
                {
                    this->statePosteriors[s][k][i][t] = alphas_2[t][k] * betas_2[t][k] / sum;
                    if (std::isnan(this->statePosteriors[s][k][i][t])) std::cout << "ERROR: statePosterior is nan! (alpha_2: " << alphas_2[t][k] << " betas_2[t][k]: " << betas_2[t][k] << " sum: " << sum << ")" << std::endl;
                }
            }

            // update init probs
            for (unsigned k = 0; k < this->K; ++k)
                this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];   

            // compute new transitioon probs
            String<String<long double> > p_i;
            resize(p_i, this->K, Exact());
            long double p_2_2_i = 0.0;
            long double p_2_3_i = 0.0;
            for (unsigned k_1 = 0; k_1 < this->K; ++k_1)    
            {
                resize(p_i[k_1], this->K, Exact());
                for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                {
                    p_i[k_1][k_2] = 0.0;
                    for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
                    {
                        p_i[k_1][k_2] += alphas_2[t-1][k_1] * this->transMatrix[k_1][k_2] * this->eProbs[s][i][t][k_2] * betas_2[t][k_2];
                        // learn p[2->2/3] for region over nThresholdForP
                        if (k_1 == 2 && k_2 == 2 && setObs[s][i].nEstimates[t] >= options.nThresholdForTransP)
                            p_2_2_i += alphas_2[t-1][k_1] * this->transMatrix[k_1][k_2] * this->eProbs[s][i][t][k_2] * betas_2[t][k_2];

                        if (k_1 == 2 && k_2 == 3 && setObs[s][i].nEstimates[t] >= options.nThresholdForTransP)
                            p_2_3_i += alphas_2[t-1][k_1] * this->transMatrix[k_1][k_2] * this->eProbs[s][i][t][k_2] * betas_2[t][k_2];
                    }
                    SEQAN_OMP_PRAGMA(critical)
                    p[k_1][k_2] += p_i[k_1][k_2];
                }
            }
            SEQAN_OMP_PRAGMA(critical)
            p_2_2 += p_2_2_i;
            SEQAN_OMP_PRAGMA(critical)
            p_2_3 += p_2_3_i;
        }
        if (stop) return false;
    }
    // update transition matrix
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
    {
        long double denumerator = 0.0;
        for (unsigned k_3 = 0; k_3 < this->K; ++k_3)
        {
            denumerator += p[k_1][k_3]; 
        }

        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
        {
            A[k_1][k_2] = p[k_1][k_2] / denumerator;

            if (A[k_1][k_2] <= 0.0) 
                A[k_1][k_2] = DBL_MIN;          // make sure not getting zero
        }
    }
    // Fix p[2->2/3] using only transition prob for region over nThresholdForP, while keeping sum of p[2->2] and p[2->3] constant! (if user options says so) 
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

        if (A[3][3] < options.minTransProbCS)
            A[3][3] = options.minTransProbCS;
        std::cout << "NOTE: Prevented transition probability '2' -> '3' from dropping below min. value of " << options.minTransProbCS << ". Set for transitions '2' -> '3' (and if necessary also for '3'->'3') to " << options.minTransProbCS << "." << std::endl;
    }
    this->transMatrix = A;
    return true;
}

// without updating transition probabilities 
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::computeStatePosteriorsFB(AppOptions &options)
{
    String<String<long double> > A = this->transMatrix;
    String<String<long double> > p;
    resize(p, this->K, Exact());
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)    
    {
        SEQAN_OMP_PRAGMA(critical)
        resize(p[k_1], this->K, Exact());
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            p[k_1][k_2] = 0.0;
    }

    for (unsigned s = 0; s < 2; ++s)
    {
        bool stop = false;
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            if (!this->setObs[s][i].discard)
            {
                unsigned T = setObs[s][i].length();
                // forward probabilities
                String<String<TDOUBLE> > alphas_1;
                String<String<TDOUBLE> > alphas_2;
                resize(alphas_1, T, Exact());
                resize(alphas_2, T, Exact());
                for (unsigned t = 0; t < T; ++t)
                {
                    resize(alphas_1[t], this->K, Exact());
                    resize(alphas_2[t], this->K, Exact());
                } 
                if (!iForward(alphas_1, alphas_2, s, i, options))
                {
                    stop = true;
                    continue;
                }
                // backward probabilities  
                String<String<TDOUBLE> > betas_2;
                resize(betas_2, T, Exact());
                for (unsigned t = 0; t < T; ++t)
                    resize(betas_2[t], this->K, Exact());
                if (!iBackward(betas_2, alphas_1, s, i))
                {
                    stop = true;
                    continue;
                }

                // compute state posterior probabilities
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                {
                    long double sum = 0.0;
                    for (unsigned k = 0; k < this->K; ++k)
                        sum += alphas_2[t][k] * betas_2[t][k];

                    if (sum == 0.0) 
                    {
                        std::cout << "ERROR: while computing state posterior probabilities! Sum of alpha and beta values is 0 at i: " << i << " t: "<< t << "." << std::endl;
                        if (options.verbosity >= 2)
                        {
                            for (unsigned k = 0; k < this->K; ++k)
                            {
                                std::cout << "k: " << k << std::endl;
                                std::cout << "alphas_2[k]: " << alphas_2[t][k] << " betas_2[t][k]: " << betas_2[t][k] << std::endl;
                            }
                        }
                        stop = true;
                        continue;
                    }
                    for (unsigned k = 0; k < this->K; ++k)
                        this->statePosteriors[s][k][i][t] = alphas_2[t][k] * betas_2[t][k] / sum;
                }

                // update init probs
                for (unsigned k = 0; k < this->K; ++k)
                    this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];   
            }
        }
        if (stop) return false;
    }
    return true;
}



template<typename TDOUBLE>
bool updateDensityParams2(String<String<String<TDOUBLE> > > &statePosteriors1, String<String<String<TDOUBLE> > > &statePosteriors2, String<String<Observations> > &setObs, GAMMA2<TDOUBLE> &d1, GAMMA2<TDOUBLE> &d2, AppOptions &options)
{
    if (!d1.updateThetaAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
        return false;

    if (!d2.updateThetaAndK(statePosteriors2, setObs, options.g2_kMin, options.g2_kMax, options))         // make sure g1k <= g2k
        return false;

    // make sure gamma1.mu < gamma2.mu   
    checkOrderG1G2(d1, d2, options);
    return true;
}

template<typename TDOUBLE>
bool updateDensityParams2(String<String<String<TDOUBLE> > > &statePosteriors1, String<String<String<TDOUBLE> > > &statePosteriors2, String<String<Observations> > &setObs, GAMMA2_REG<TDOUBLE> &d1, GAMMA2_REG<TDOUBLE> &d2, AppOptions &options)
{
    if (!d1.updateRegCoeffsAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
        return false;

    double g2_kMin = options.g2_kMin;
    if (options.g1_k_le_g2_k)
        g2_kMin = std::max(d1.k, options.g2_kMin);

    if (!d2.updateRegCoeffsAndK(statePosteriors2, setObs, g2_kMin, options.g2_kMax, options))
        return false;

    // make sure gamma1.mu < gamma2.mu    
    checkOrderG1G2(d1, d2, options);
    return true;
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::updateDensityParams(TGAMMA &d1, TGAMMA &d2, AppOptions &options)   
{
    String<String<String<TDOUBLE> > > statePosteriors1;
    String<String<String<TDOUBLE> > > statePosteriors2;
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

    updateDensityParams2(statePosteriors1, statePosteriors2, this->setObs, d1, d2, options); 

    return true;
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool HMM<TGAMMA, TBIN, TDOUBLE>::updateDensityParams(TGAMMA /*&d1*/, TGAMMA /*&d2*/, TBIN &bin1, TBIN &bin2, AppOptions &options)   
{
    String<String<String<TDOUBLE> > > statePosteriors1;
    String<String<String<TDOUBLE> > > statePosteriors2;
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
    bin1.updateP(statePosteriors1, this->setObs, options); 
    bin2.updateP(statePosteriors2, this->setObs, options);

    // make sure bin1.p < bin2.p   
    checkOrderBin1Bin2(bin1, bin2);

    return true;
}



// Baum-Welch
// with scaling
template<typename TGAMMA, typename TBIN, typename TDOUBLE> 
bool HMM<TGAMMA, TBIN, TDOUBLE>::baumWelch(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, CharString learnTag, AppOptions &options)
{
    TGAMMA prev_d1 = d1;
    TGAMMA prev_d2 = d2;
    TBIN prev_bin1 = bin1;
    TBIN prev_bin2 = bin2;
    for (unsigned iter = 0; iter < options.maxIter_bw; ++iter)
    {
        std::cout << ".. " << iter << "th iteration " << std::endl;
        std::cout << "                        computeEmissionProbs() " << std::endl;
        if (!computeEmissionProbs(d1, d2, bin1, bin2, true, options) )
        {
            std::cerr << "ERROR: Could not compute emission probabilities! " << std::endl;
            return false;
        }
        std::cout << "                        computeStatePosteriorsFB() " << std::endl;
        if (!computeStatePosteriorsFBupdateTrans(options))
            return false;
        
        std::cout << "                        updateDensityParams() " << std::endl;

        if (learnTag == "LEARN_BINOMIAL")
        {
            if (!updateDensityParams(d1, d2, bin1, bin2, options))
            {
                std::cerr << "ERROR: Could not update parameters! " << std::endl;
                return false;
            }
        }
        else
        {
            if (!updateDensityParams(d1, d2, options))
            {
                std::cerr << "ERROR: Could not update parameters! " << std::endl;
                return false;
            }
        }
        
        if (learnTag == "LEARN_GAMMA" && checkConvergence(d1, prev_d1, options) && checkConvergence(d2, prev_d2, options) )             
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        else if (learnTag != "LEARN_GAMMA" && checkConvergence(bin1, prev_bin1, options) && checkConvergence(bin2, prev_bin2, options) )             
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        prev_d1 = d1;
        prev_d2 = d2;
        prev_bin1 = bin1;
        prev_bin2 = bin2;

        myPrint(d1);
        myPrint(d2);
        
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
            myPrint(bin1);
            myPrint(bin2);
        }
    }
    return true;
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE> 
bool HMM<TGAMMA, TBIN, TDOUBLE>::applyParameters(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, AppOptions &options)
{
    if (!computeEmissionProbs(d1, d2, bin1, bin2, false, options))
    {
        std::cerr << "ERROR: Could not compute emission probabilities! " << std::endl;
        return false;
    }
    if (!computeStatePosteriorsFB(options))
        return false;

    return true;
}


// returns log P
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
long double HMM<TGAMMA, TBIN, TDOUBLE>::viterbi(String<String<String<__uint8> > > &states)
{
    long double p = 1.0;
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(states[s], length(this->setObs[s]), Exact());
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            if (!this->setObs[s][i].discard)
            {
                resize(states[s][i], this->setObs[s][i].length(), Exact());
                // store for each t and state maximizing precursor joint probability of state sequence and observation
                String<String<long double> > vits;
                resize(vits, this->setObs[s][i].length(), Exact());
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                    resize(vits[t], this->K, Exact());
                // store for each t and state maximizing precursor state
                String<String<unsigned> > track;
                resize(track, this->setObs[s][i].length(), Exact());
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                    resize(track[t], this->K, Exact());

                // initialize
                for (unsigned k = 0; k < this->K; ++k)
                    vits[0][k] = this->initProbs[s][i][k] * this->eProbs[s][i][0][k];
                // recursion
                for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
                {
                    for (unsigned k = 0; k < this->K; ++k)
                    {
                        long double max_v = vits[t-1][0] * this->transMatrix[0][k];
                        unsigned max_k = 0;
                        for (unsigned k_p = 1; k_p < this->K; ++k_p)
                        {
                            long double v = vits[t-1][k_p] * this->transMatrix[k_p][k];
                            if (v > max_v)
                            {
                                max_v = v;
                                max_k = k_p;
                            }
                        }
                        vits[t][k] = max_v * this->eProbs[s][i][t][k];
                        track[t][k] = max_k;
                    }
                }
                // backtracking
                long double max_v = vits[this->setObs[s][i].length() - 1][0];
                unsigned max_k = 0;
                for (unsigned k = 1; k < this->K; ++k)
                {
                    if (vits[this->setObs[s][i].length() - 1][k] >= max_v)
                    {
                        max_v = vits[this->setObs[s][i].length() - 1][k];
                        max_k = k;
                    }
                }
                states[s][i][this->setObs[s][i].length() - 1] = max_k;
                for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
                    states[s][i][t] = track[t+1][states[s][t+1]];

                p *= max_v;
            }
        }
    }
    return p;   
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
long double HMM<TGAMMA, TBIN, TDOUBLE>::viterbi_log(String<String<String<__uint8> > > &states)
{
    long double p = 0.0;
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(states[s], length(this->setObs[s]), Exact());
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            if (!this->setObs[s][i].discard)
            {
                resize(states[s][i], this->setObs[s][i].length(), Exact());
                // store for each t and state maximizing precursor joint probability of state sequence and observation
                String<String<long double> > vits;
                resize(vits, this->setObs[s][i].length(), Exact());
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                    resize(vits[t], this->K, Exact());
                // store for each t and state maximizing precursor state
                String<String<unsigned> > track;
                resize(track, this->setObs[s][i].length(), Exact());
                for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
                    resize(track[t], this->K, Exact());

                // SEQAN_ASSERT_GT( ,0.0) or <- DBL_MIN

                // initialize
                for (unsigned k = 0; k < this->K; ++k)
                    vits[0][k] = log(this->initProbs[s][i][k]) + log(this->eProbs[s][i][0][k]);
                // recursion
                for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
                {
                    for (unsigned k = 0; k < this->K; ++k)
                    {
                        long double max_v = vits[t-1][0] + log(this->transMatrix[0][k]);
                        unsigned max_k = 0;
                        for (unsigned k_p = 1; k_p < this->K; ++k_p)
                        {
                            double v = vits[t-1][k_p] + log(this->transMatrix[k_p][k]);
                            if (v > max_v)
                            {
                                max_v = v;
                                max_k = k_p;
                            }
                        }
                        vits[t][k] = max_v + log(this->eProbs[s][i][t][k]);
                        track[t][k] = max_k;
                    }
                }
                // backtracking
                long double max_v = vits[this->setObs[s][i].length() - 1][0];
                unsigned max_k = 0;
                for (unsigned k = 1; k < this->K; ++k)
                {
                    if (vits[this->setObs[s][i].length() - 1][k] >= max_v)
                    {
                        max_v = vits[this->setObs[s][i].length() - 1][k];
                        max_k = k;
                    }
                }
                states[s][i][this->setObs[s][i].length() - 1] = max_k;
                for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
                    states[s][i][t] = track[t+1][states[s][i][t+1]];

                p += max_v;
            }
        }
    }
    // NOTE p: of all sites, not only selected for parameter fitting, not necessarly increases!
    return p;        
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
void HMM<TGAMMA, TBIN, TDOUBLE>::posteriorDecoding(String<String<String<__uint8> > > &states)
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
                    if (s==0 && (t + this->setPos[s][i] - 1 == 22418895) )
                    {
                        std::cout << "State posterior probs" << std::endl;
                        std::cout << "k: 0  " << this->statePosteriors[s][0][i][t] << std::endl;
                        std::cout << "k: 1  " << this->statePosteriors[s][1][i][t] << std::endl;
                        std::cout << "k: 2  " << this->statePosteriors[s][2][i][t] << std::endl;
                        std::cout << "k: 3  " << this->statePosteriors[s][3][i][t] << std::endl;
                        std::cout << "max k  " << (int)states[s][i][t] << std::endl;
                        
                    }
                }
            }
        }
    }
}


template<typename TDOUBLE>
void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, String<String<Observations> > &setObs, GAMMA2_REG<TDOUBLE> &g1)
{
    double b0 = g1.b0;
    double b1 = g1.b1;
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(setObs[s]); ++i)
        {
            if (!setObs[s][i].discard)
            {
                for (unsigned t = 0; t < setObs[s][i].length(); ++t)
                {
                    double x1 = setObs[s][i].rpkms[t];
                    double g1_pred = exp(b0 + b1 * x1);
                    if (states[s][i][t] >= 2 && setObs[s][i].kdes[t] < g1_pred)
                        states[s][i][t] -= 2;
                }
            }
        }
    }
}

template<typename TDOUBLE>
void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, String<String<Observations> > &setObs, GAMMA2<TDOUBLE> &g1)
{
    // do nothing
}

// for GLM with input signal: 
// when using free gamma shapes, i.e. gamma1.k can be > gamma2.k
// make sure sites with fragment coverage (KDE) below gamma1.mean are classified as 'non-enriched'
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
void HMM<TGAMMA, TBIN, TDOUBLE>::rmBoarderArtifacts(String<String<String<__uint8> > > &states, TGAMMA &g1)
{
    rmBoarderArtifacts2(states, this->setObs, g1);
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
                        ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, min_val) );
                    }
                    else if (options.score_type == 1)
                    {
                        // log(3/2)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/data.statePosteriors[s][2][i][t]);
                    }
                    else if (options.score_type == 2)
                    {
                        // log(3/1)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/data.statePosteriors[s][1][i][t]);
                    }
                    else if (options.score_type == 3)
                    {
                        // log(3/0)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/data.statePosteriors[s][0][i][t]);
                    }                    
                    else if (options.score_type == 4)
                    {
                        // log(3/(1+2))  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/(data.statePosteriors[s][1][i][t] + data.statePosteriors[s][2][i][t]));
                    }
                    else if (options.score_type == 5)
                    {
                        // 3 
                        ss << (double)data.statePosteriors[s][3][i][t];
                    }
                    else if (options.score_type == 6)
                    {
                        // log(enriched/non-enriched) + log(crosslinked/non-crosslinked)  
                        ss << (double)log((data.statePosteriors[s][2][i][t] + data.statePosteriors[s][3][i][t])/(data.statePosteriors[s][0][i][t] + data.statePosteriors[s][1][i][t])) + (double)log((data.statePosteriors[s][1][i][t] + data.statePosteriors[s][3][i][t])/(data.statePosteriors[s][0][i][t] + data.statePosteriors[s][2][i][t])) ;
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
                        ss << (double)log(data.statePosteriors[s][3][i][t]/data.statePosteriors[s][2][i][t]);
                    }
                    else if (options.score_type == 2)
                    {
                        // log(3/1)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/data.statePosteriors[s][1][i][t]);
                    }
                    else if (options.score_type == 3)
                    {
                        // log(3/0)  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/data.statePosteriors[s][0][i][t]);
                    }                    
                    else if (options.score_type == 4)
                    {
                        // log(3/(1+2))  
                        ss << (double)log(data.statePosteriors[s][3][i][t]/(data.statePosteriors[s][1][i][t] + data.statePosteriors[s][2][i][t]));
                    }
                    else if (options.score_type == 5)
                    {
                        // 3 
                        ss << (double)data.statePosteriors[s][3][i][t];
                    }
                    else if (options.score_type == 6)
                    {
                        // log(enriched/non-enriched) + log(crosslinked/non-crosslinked)  
                        ss << (double)log((data.statePosteriors[s][2][i][t] + data.statePosteriors[s][3][i][t])/(data.statePosteriors[s][0][i][t] + data.statePosteriors[s][1][i][t])) + (double)log((data.statePosteriors[s][1][i][t] + data.statePosteriors[s][3][i][t])/(data.statePosteriors[s][0][i][t] + data.statePosteriors[s][2][i][t])) ;
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

  
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
void myPrint(HMM<TGAMMA, TBIN, TDOUBLE> &hmm)
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
