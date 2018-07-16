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

    __uint8                    K;                  // no. of sates
    String<String<String<double> > >   initProbs;          // intital probabilities

    String<String<Observations> >       & setObs;          // workaround for partial specialization
    String<String<unsigned> >           & setPos;
    unsigned                            contigLength;
    String<String<double> >     transMatrix;

    HMM(int K_, String<String<Observations> > & setObs_, String<String<unsigned> > & setPos_, unsigned &contigLength_): K(K_), setObs(setObs_), setPos(setPos_), contigLength(contigLength)
    {
        // initialize transition probabilities
        resize(transMatrix, K, Exact());
        double trans1 = 0.6;    // increa:sed probability to stay in same state
        for (unsigned i = 0; i < K; ++i)
        {
            resize(transMatrix[i], K, Exact());
            for (unsigned j = 0; j < K; ++j)
                if (i == j)
                    transMatrix[i][j] = trans1;
                else
                    transMatrix[i][j] = (1.0 - trans1) / (K - 1.0);
        }
       
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
                    initProbs[s][i][k] = 1.0 / K;

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

    HMM<TGAMMA, TBIN, TDOUBLE>();
    ~HMM<TGAMMA, TBIN, TDOUBLE>();
    void setInitProbs(String<double> &probs);
    bool computeEmissionProbs(TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2, bool learning, AppOptions &options);
    bool iForward(String<String<TDOUBLE> > &alphas_1, String<String<TDOUBLE> > &alphas_2, unsigned s, unsigned i, AppOptions &options);
    //void forward_noSc();
    bool iBackward(String<String<TDOUBLE> > &betas_2, String<String<TDOUBLE> > &alphas_1, unsigned s, unsigned i);
    //void backward_noSc();
    bool computeStatePosteriorsFB(AppOptions &options);
    bool computeStatePosteriorsFBupdateTrans(AppOptions &options);
    //void updateTransition(AppOptions &options);
    //void updateTransition2();
    //void updateTransition_noSc2();
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
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]);
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    // 
    if (eProbs[0] == 0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0)
    {
        if (options.verbosity >= 2)
        {
            std::cout << "WARNING: all emission probabilities are 0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << d1.getDensity(setObs.kdes[t]) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << d2.getDensity(setObs.kdes[t]) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]) << std::endl;
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
    double x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
    double d1_pred = exp(d1.b0 + d1.b1 * x);
    double d2_pred = exp(d2.b0 + d2.b1 * x);

    long double g1_d = 1.0;
    long double g2_d = 0.0;
    if (setObs.kdes[t] >= d1.tp) 
    {
        g1_d = d1.getDensity(setObs.kdes[t], d1_pred);
        g2_d = d2.getDensity(setObs.kdes[t], d2_pred); 
    }
    long double  bin1_d = 1.0;
    long double  bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]);
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    // 
    if ((eProbs[0] == 0.0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0) || 
            (std::isnan(eProbs[0]) || std::isnan(eProbs[1]) || std::isnan(eProbs[2]) || std::isnan(eProbs[3])) ||
            (std::isinf(eProbs[0]) || std::isinf(eProbs[1]) || std::isinf(eProbs[2]) || std::isinf(eProbs[3])))
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "WARNING: all emission probabilities are 0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate b: " << x << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << d1.getDensity(setObs.kdes[t], d1_pred) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << d2.getDensity(setObs.kdes[t], d2_pred) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t]) << std::endl;
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
    double bin1_pred = 1.0/(1.0+exp(-bin1.b0 - bin1.regCoeffs[mId]*setObs.fimoScores[t]));
    double bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*setObs.fimoScores[t]));

    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred);
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    //
    if (eProbs[0] == 0.0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0)
    {
        if (options.verbosity >= 2)
        {
            std::cout << "WARNING: all emission probabilities are 0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate x: " << setObs.fimoScores[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << d1.getDensity(setObs.kdes[t]) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << d2.getDensity(setObs.kdes[t]) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred) << std::endl;
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
    double x = std::max(setObs.rpkms[t], options.minRPKMtoFit);
    double d1_pred = exp(d1.b0 + d1.b1 * x);
    double d2_pred = exp(d2.b0 + d2.b1 * x);

    long double g1_d = 1.0;
    long double g2_d = 0.0;
    if (setObs.kdes[t] >= d1.tp) 
    {
        g1_d = d1.getDensity(setObs.kdes[t], d1_pred);
        g2_d = d2.getDensity(setObs.kdes[t], d2_pred); 
    }
    unsigned mId = setObs.motifIds[t];
    double bin1_pred = 1.0/(1.0+exp(-bin1.b0 - bin1.regCoeffs[mId]*setObs.fimoScores[t]));
    double bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*setObs.fimoScores[t]));

    long double bin1_d = 1.0;
    long double bin2_d = 0.0;
    if (setObs.truncCounts[t] > 0)
    {
        bin1_d = bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred);
        bin2_d = bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred);
    }
    eProbs[0] = g1_d * bin1_d;    
    eProbs[1] = g1_d * bin2_d;
    eProbs[2] = g2_d * bin1_d;
    eProbs[3] = g2_d * bin2_d;

    // 
    if ((eProbs[0] == 0.0 && eProbs[1] == 0.0 && eProbs[2] == 0.0 && eProbs[3] == 0.0) || 
            (std::isnan(eProbs[0]) || std::isnan(eProbs[1]) || std::isnan(eProbs[2]) || std::isnan(eProbs[3])) ||
            (std::isinf(eProbs[0]) || std::isinf(eProbs[1]) || std::isinf(eProbs[2]) || std::isinf(eProbs[3])))
    {
        if (options.verbosity >= 2)
        {
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "WARNING: all emission probabilities are 0!" << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       fragment coverage (kde): " << setObs.kdes[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       read start count: " << (int)setObs.truncCounts[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       estimated n: " << setObs.nEstimates[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate b: " << x << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       covariate x: " << setObs.fimoScores[t] << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-enriched' gamma: " << d1.getDensity(setObs.kdes[t], d1_pred) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'enriched' gamma: " << d2.getDensity(setObs.kdes[t], d2_pred) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'non-crosslink' binomial: " << bin1.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin1_pred) << std::endl;
            SEQAN_OMP_PRAGMA(critical) 
                std::cout << "       emission probability 'crosslink' binomial: " << bin2.getDensity(setObs.truncCounts[t], setObs.nEstimates[t], bin2_pred) << std::endl;
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
                std::cout << "ERROR: Emission probability became 0.0! This might be due to a too small learning subset or due to some extreme artifacts." << std::endl;
                SEQAN_OMP_PRAGMA(critical) 
                if (s == 0) 
                    std::cout << " Interval: [" << (this->setPos[s][i]) << ", " << (this->setPos[s][i] + this->setObs[s][i].length()) << ") on forward strand." << std::endl;
                else 
                    std::cout << " Interval: [" << (this->contigLength - this->setPos[s][i] - 1) << ", " << (this->contigLength - this->setPos[s][i] - 1 + this->setObs[s][i].length()) << ") on reverse strand." << std::endl;
                stop = true;
                if (!options.useHighPrecision)
                {
                    SEQAN_OMP_PRAGMA(critical) 
                    std::cout << "NOTE: Try running PureCLIP with parameter '-ld' ()!" << std::endl;
                }
            }
            else if (!learning && discardInterval)
            {
                SEQAN_OMP_PRAGMA(critical) 
                if (s == 0) 
                    std::cout << "Warning: discarding interval [" << (this->setPos[s][i]) << ", " << (this->setPos[s][i] + this->setObs[s][i].length()) << ") on forward strand due to emission probabilities of 0 (set to state 'non-enriched + non-crosslink')." << std::endl;
                else 
                    std::cout << "Warning: discarding interval [" << (this->contigLength - this->setPos[s][i] - 1) << ", " << (this->contigLength - this->setPos[s][i] - 1 + this->setObs[s][i].length()) << ") on reverse strand due to emission probabilities of 0 (set to state 'non-enriched + non-crosslink')." << std::endl;
                if (!options.useHighPrecision)
                {
                    SEQAN_OMP_PRAGMA(critical) 
                    std::cout << "NOTE: If this happens frequently, rerun PureCLIP with parameter '-ld' ()!" << std::endl;
                }
            }
        }   
    }
    if (stop) return false;
    return true;
}


 
// Forward-backward algorithm 
// without scaling for testing 
/*template<typename TD1, typename TD2, typename TB1, typename TB2>
void HMM<TGAMMA, TBIN, TDOUBLE>::forward_noSc()
{
    for (unsigned s = 0; s < 2; ++s)
    {
        // for t = 1
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            for (unsigned k = 0; k < this->K; ++k)
                this->alphas_2[s][i][0][k] = this->initProbs[s][i][k] * this->eProbs[s][i][0][k];
        
            // for t = 2:T
            for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
            {
                for (unsigned k = 0; k < this->K; ++k)
                {
                    // sum over previous states
                    double sum = 0.0;
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        sum += alphas_2[s][i][t-1][k_2] * this->transMatrix[k_2][k];

                    this->alphas_2[s][i][t][k] = sum * this->eProbs[s][i][t][k];
                }
            }
        }
    }
}

//// probabilities getting small, in a naive implementation alpha and beta will underflow
//// -> scaling technique (more stable and faster than log method ?)
// compute scaled alphas for each state and corresponding c values
template<typename TD1, typename TD2, typename TB1, typename TB2>
void HMM<TGAMMA, TBIN, TDOUBLE>::forward()
{
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            // for t = 1
            double norm = 0.0;
            for (unsigned k = 0; k < this->K; ++k)
            {
                this->alphas_1[s][i][0][k] = this->initProbs[s][i][k] * this->eProbs[s][i][0][k];
                norm += this->alphas_1[s][i][0][k];
            }
            if (norm == 0.0) 
            {
                std::cerr << "ERROR: norm == 0 at t: "<< 0 << "  i: " << i << std::endl;
                for (unsigned k = 0; k < this->K; ++k)
                {
                    std::cout << "k: " << k << std::endl;
                    std::cout << "this->eProbs[s][i][0][k] " << this->eProbs[s][i][0][k] << std::endl;
                }
            }

            for (unsigned k = 0; k < this->K; ++k)
               this->alphas_2[s][i][0][k] = alphas_1[s][i][0][k] / norm;

            // for t = 2:T
            for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
            {
                norm = 0.0;
                for (unsigned k = 0; k < this->K; ++k)
                {
                    // sum over previous states
                    double sum = 0.0;
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        sum += alphas_2[s][i][t-1][k_2] * this->transMatrix[k_2][k];
                    
                    if (sum == 0.0 || std::isnan(sum)) 
                    {
                        std::cerr << "ERROR: sum = " << norm << " at t: "<< t << "  i: " << i << std::endl;
                        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                            std::cout << " k_2: " << k_2 <<  " alphas_2[s][i][t-1][k_2]: " << alphas_2[s][i][t-1][k_2] << " transMatrix[k_2][k]: " << this->transMatrix[k_2][k] << std::endl; 
                    }

                    // alpha_1
                    this->alphas_1[s][i][t][k] = sum * this->eProbs[s][i][t][k];        // - nan?
                    norm += this->alphas_1[s][i][t][k];
                }
                
                if (norm == 0.0 || std::isnan(norm)) 
                {
                    std::cerr << "ERROR: norm = " << norm << " at t: "<< t << "  i: " << i << std::endl;
                    for (unsigned k = 0; k < this->K; ++k)
                    {
                        std::cout << "k: " << k << std::endl;
                        std::cout << "this->eProbs[s][i][t][k] " << this->eProbs[s][i][t][k] << std::endl;
                    }
                }
                // normalize
                for (unsigned k = 0; k < this->K; ++k)
                    this->alphas_2[s][i][t][k] = alphas_1[s][i][t][k] / norm;   // TODO store scaling coefficients too !?
            }
        }
    }
}*/


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
// without scaling for testing 
/*template<typename TD1, typename TD2, typename TB1, typename TB2>
void HMM<TGAMMA, TBIN, TDOUBLE>::backward_noSc()
{
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            // for t = T
            for (unsigned k = 0; k < this->K; ++k)
               this->betas_2[s][i][this->setObs[s][i].length() - 1][k] = 1.0;
            
            // for t = 2:T
            for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
            {
                for (unsigned k = 0; k < this->K; ++k)
                {
                    // sum over previous states
                    double sum = 0.0;
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        sum += betas_2[s][i][t+1][k_2] * this->transMatrix[k][k_2] * this->eProbs[s][i][t+1][k_2];
                    // beta_1
                    this->betas_2[s][i][t][k] = sum;
                }
            }
        }
    }
}

// with scaling method
template<typename TD1, typename TD2, typename TB1, typename TB2>
void HMM<TGAMMA, TBIN, TDOUBLE>::backward()
{
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif  
        for (unsigned i = 0; i < length(this->setObs[s]); ++i)
        {
            // for t = T
            for (unsigned k = 0; k < this->K; ++k)
               this->betas_1[s][i][this->setObs[s][i].length() - 1][k] = 1.0;
            
            double norm = 0.0;      // use scaling coefficients from alphas here !
            for (unsigned k = 0; k < this->K; ++k)
               norm += this->alphas_1[s][i][this->setObs[s][i].length() - 1][k];

            for (unsigned k = 0; k < this->K; ++k)
               this->betas_2[s][i][this->setObs[s][i].length() - 1][k] = betas_1[s][i][this->setObs[s][i].length() - 1][k] / norm;

            // for t = 2:T
            for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
            {
                norm = 0.0;
                for (unsigned k = 0; k < this->K; ++k)      // precompute ???
                    norm += this->alphas_1[s][i][t][k];

                for (unsigned k = 0; k < this->K; ++k)
                {
                    // sum over previous states
                    double sum = 0.0;
                    for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                        sum += betas_2[s][i][t+1][k_2] * this->transMatrix[k][k_2] * this->eProbs[s][i][t+1][k_2];
                    
                    // beta_1
                    this->betas_1[s][i][t][k] = sum;
                    // beta_2
                    this->betas_2[s][i][t][k] = this->betas_1[s][i][t][k] / norm;
                }
            }
        }
    }
}*/


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

    for (unsigned k = 0; k < this->K; ++k)
       betas_2[this->setObs[s][i].length() - 1][k] = betas_1[this->setObs[s][i].length() - 1][k] / norm;

    // for t = 2:T
    for (int t = this->setObs[s][i].length() - 2; t >= 0; --t)
    {
        norm = 0.0;
        for (unsigned k = 0; k < this->K; ++k)      // precompute ???
            norm += alphas_1[t][k];

        for (unsigned k = 0; k < this->K; ++k)
        {
            // sum over previous states
            long double sum = 0.0;
            for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                sum += betas_2[t+1][k_2] * this->transMatrix[k][k_2] * this->eProbs[s][i][t+1][k_2];
            
            // beta_1
            betas_1[t][k] = sum;
            // beta_2
            betas_2[t][k] = betas_1[t][k] / norm;
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
    String<String<double> > A = this->transMatrix;
    String<String<double> > p;
    resize(p, this->K, Exact());
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)    
    {
        SEQAN_OMP_PRAGMA(critical)
        resize(p[k_1], this->K, Exact());
        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
            p[k_1][k_2] = 0.0;
    }
    double p_2_2 = 0.0;     // for separate learning of trans. prob from '2' -> '2'
    double p_2_3 = 0.0;     // for separate learning of trans. prob from '2' -> '3' 

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
                    this->statePosteriors[s][k][i][t] = alphas_2[t][k] * betas_2[t][k] / sum;
            }

            // update init probs
            for (unsigned k = 0; k < this->K; ++k)
                this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];   

            // compute new transitioon probs
            String<String<double> > p_i;
            resize(p_i, this->K, Exact());
            double p_2_2_i = 0.0;
            double p_2_3_i = 0.0;
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
        double denumerator = 0.0;
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
        double sum_2_23 = A[2][2] + A[2][3];
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
    String<String<double> > A = this->transMatrix;
    String<String<double> > p;
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
        if (stop) return false;
    }
    return true;
}


// both for scaling and no-scaling
/*template<typename TD1, typename TD2, typename TB1, typename TB2>
void HMM<TGAMMA, TBIN, TDOUBLE>::updateTransition(AppOptions &options)    // precompute numerator, denumerator, for each k_1, k_2 combination!!!
{
    String<String<double> > A = this->transMatrix;
    String<String<double> > p;
    resize(p, this->K, Exact());

#if HMM_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1)) 
#endif
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)    
    {
        SEQAN_OMP_PRAGMA(critical)
        resize(p[k_1], this->K, Exact());

        for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
        {
            p[k_1][k_2] = 0.0;
            for (unsigned s = 0; s < 2; ++s)
                for (unsigned i = 0; i < length(this->setObs[s]); ++i)
                    for (unsigned t = 1; t < this->setObs[s][i].length(); ++t)
                        p[k_1][k_2] += this->alphas_2[s][i][t-1][k_1] * this->transMatrix[k_1][k_2] *  this->eProbs[s][i][t][k_2] * betas_2[s][i][t][k_2]; 
        }
    }
    
    for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
    {
        double denumerator = 0.0;
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
    // keep transProb of '2' -> '3' on min. value
    if (A[2][3] < options.minTransProbCS)
    {
        A[2][3] = options.minTransProbCS;

        if (A[3][3] < options.minTransProbCS)
            A[3][3] = options.minTransProbCS;
        // TODO normalize again ?
        std::cout << "NOTE: Prevented transition probability '2' -> '3' from dropping below min. value of " << options.minTransProbCS << ". Set for transitions '2' -> '3' (and if necessary also for '3'->'3') to " << options.minTransProbCS << "." << std::endl;
    }
    this->transMatrix = A;
}*/

// no scaling, version 2 (does this work?)
/*template<typename TD1, typename TD2>
void HMM<TD1, TD2>::updateTransition_noSc2(String<String<double> > &A, unsigned k_1, unsigned k_2)
{
    double numerator = 0.0;
    for (unsigned t = 1; t < this->T; ++t)
    {
        double sum1 = this->alphas_2[t-1][k_1] * this->transMatrix[k_1][k_2] *  this->eProbs[t][k_2] * betas_2[t][k_2];  
        double sum2 = 0.0;
        for (unsigned k = 0; k <this->K; ++k)
            sum2 += this->alphas_2[t-1][k] * betas_2[t-1][k]; 
        numerator += sum1 / sum2;
    }
    double denumerator = 0.0;
    for (unsigned t = 0; t < this->T - 1; ++t)
        denumerator += this->statePosteriors[k_2][t];
    
    A[k_1][k_2] = numerator / denumerator;  
}
// with scaling, version 2 (does this work ?)
template<typename TD1, typename TD2>
void HMM<TD1, TD2>::updateTransition2(String<String<double> > &A, unsigned k_1, unsigned k_2)
{
    double numerator = 0.0;
    for (unsigned t = 1; t < this->T; ++t)
        numerator += this->alphas_2[t-1][k_1] * this->transMatrix[k_1][k_2] *  this->eProbs[t][k_2] * betas_2[t][k_2];  
   
    double denumerator  = 0.0;
    for (unsigned t = 1; t < this->T; ++t)
    {
        double norm = 0.0;
        for (unsigned k = 0; k < this->K; ++k)      
            norm += this->alphas_1[t-1][k];

        denumerator  += this->alphas_2[t-1][k_1] * this->betas_2[t-1][k_1] * norm;
    }
    A[k_1][k_2] = numerator / denumerator ;  
}*/

template<typename TDOUBLE>
bool updateDensityParams2(String<String<String<TDOUBLE> > > &statePosteriors1, String<String<String<TDOUBLE> > > &statePosteriors2, String<String<Observations> > &setObs, GAMMA2<TDOUBLE> &d1, GAMMA2<TDOUBLE> &d2, AppOptions &options)
{
    if (options.gslSimplex2)
    {
        if (!d1.updateThetaAndK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options))
            return false;

        if (!d2.updateThetaAndK(statePosteriors2, setObs, options.g2_kMin, options.g2_kMax, options))         // make sure g1k <= g2k
            return false;

        // make sure gamma1.mu < gamma2.mu   
        checkOrderG1G2(d1, d2, options);
    }
    else    // TODO get rid of this
    {
        // 
        d1.updateTheta(statePosteriors1, setObs, options);  
        d2.updateTheta(statePosteriors2, setObs, options);

        // make sure gamma1.mu < gamma2.mu    
        checkOrderG1G2(d1, d2, options);

        d1.updateK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options);  

        d2.updateK(statePosteriors2, setObs, options.g2_kMin, options.g2_kMax, options); 
    }
    return true;
}

template<typename TDOUBLE>
bool updateDensityParams2(String<String<String<TDOUBLE> > > &statePosteriors1, String<String<String<TDOUBLE> > > &statePosteriors2, String<String<Observations> > &setObs, GAMMA2_REG<TDOUBLE> &d1, GAMMA2_REG<TDOUBLE> &d2, AppOptions &options)
{
    if (options.gslSimplex2)
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
    }
    else
    { 
        d1.updateMean(statePosteriors1, setObs, options); 
        d2.updateMean(statePosteriors2, setObs, options);

        // make sure gamma1.mu < gamma2.mu    
        checkOrderG1G2(d1, d2, options);

        d1.updateK(statePosteriors1, setObs, options.g1_kMin, options.g1_kMax, options); 

        double g2_kMin = options.g2_kMin;
        if (options.g1_k_le_g2_k)
            g2_kMin = std::max(d1.k, options.g2_kMin);

        d2.updateK(statePosteriors2, setObs, g2_kMin, options.g2_kMax, options);
    }
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
// E: compute state posterior (gamma), transition posterior (xi)
// M: estimate parameters
/*template<typename TD1, typename TD2, typename TB1, typename TB2> 
bool HMM<TGAMMA, TBIN, TDOUBLE>::baumWelch_noSc(TD1 &d1, TD2 &d2, TB1 &bin1, TB2 &bin2, CharString learnTag, AppOptions &options)
{
    double prev_p = 666.0;
    for (unsigned iter = 0; iter < options.maxIter_bw; ++iter)
    {
        std::cout << ".. " << iter << "th iteration " << std::endl;
        if (!computeEmissionProbs(d1, d2, bin1, bin2, options))
        {
            std::cerr << "ERROR: Could not compute emission probabilities! " << std::endl;
            return false;
        }
        // E-step
        forward_noSc();
         //Note: likelihood of all sites, not only selected for parameter fitting, not necessarly increases! 

        backward_noSc();
        computeStatePosteriors();
        // M-step 
        for (unsigned s = 0; s < 2; ++s)
            for (unsigned i = 0; i < length(this->setObs[s]); ++i)
                for (unsigned k = 0; k < this->K; ++k)
                    this->initProbs[s][i][k] = this->statePosteriors[s][k][i][0];   
 
        updateTransition(options);

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
    }
    return true;
}*/

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
            resize(states[s][i], this->setObs[s][i].length(), Exact());
            // store for each t and state maximizing precursor joint probability of state sequence and observation
            String<String<TDOUBLE> > vits;
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
                    TDOUBLE max_v = vits[t-1][0] * this->transMatrix[0][k];
                    unsigned max_k = 0;
                    for (unsigned k_p = 1; k_p < this->K; ++k_p)
                    {
                        TDOUBLE v = vits[t-1][k_p] * this->transMatrix[k_p][k];
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
            TDOUBLE max_v = vits[this->setObs[s][i].length() - 1][0];
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
            resize(states[s][i], this->setObs[s][i].length(), Exact());
            // store for each t and state maximizing precursor joint probability of state sequence and observation
            String<String<TDOUBLE> > vits;
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
                    TDOUBLE max_v = vits[t-1][0] + log(this->transMatrix[0][k]);
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
            TDOUBLE max_v = vits[this->setObs[s][i].length() - 1][0];
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
            resize(states[s][i], this->setObs[s][i].length(), Exact());
            for (unsigned t = 0; t < this->setObs[s][i].length(); ++t)
            {
                double max_p = 0.0;
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


template<typename TDOUBLE>
void rmBoarderArtifacts2(String<String<String<__uint8> > > &states, String<String<Observations> > &setObs, GAMMA2_REG<TDOUBLE> &g1)
{
    double b0 = g1.b0;
    double b1 = g1.b1;
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(setObs[s]); ++i)
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


void writeStates(BedFileOut &outBed,
                 Data &data,
                 FragmentStore<> &store, 
                 unsigned contigId,
                 AppOptions &options)          
{  
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.states[s]); ++i)
        {
            for (unsigned t = 0; t < length(data.states[s][i]); ++t)
            {
                if (options.outputAll && data.setObs[s][i].truncCounts[t] >= 1)
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

                    // log posterior prob. ratio score
                    double secondBest = 0.0;
                    for (unsigned k = 0; k < 4; ++k)
                    {
                        if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                            secondBest = data.statePosteriors[s][k][i][t];
                    }                    
                    ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) );

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

                    writeRecord(outBed, record);
                }
                else if (data.states[s][i][t] == 3)
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

                    // log posterior prob. ratio score
                    double secondBest = 0.0;
                    for (unsigned k = 0; k < 4; ++k)
                    {
                        if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                            secondBest = data.statePosteriors[s][k][i][t];
                    }                    
                    ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) );

                    record.score = ss.str();
                    ss.str("");  
                    ss.clear();  
                    if (s == 0)
                        record.strand = '+';
                    else
                        record.strand = '-';

                    writeRecord(outBed, record);
                }               
            }
        }
    }
}


void writeRegions(BedFileOut &outBed,
                 Data &data,
                 FragmentStore<> &store, 
                 unsigned contigId,
                 AppOptions &options)          
{  
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.states[s]); ++i)
        {
            for (unsigned t = 0; t < length(data.states[s][i]); ++t)
            {
                if (data.states[s][i][t] == 3)
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
                    double secondBest = 0.0;
                    for (unsigned k = 0; k < 4; ++k)
                    {
                        if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                            secondBest = data.statePosteriors[s][k][i][t];
                    }                    
                    ss << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) );
                    record.score = ss.str();
                    ss.str("");  
                    ss.clear();  
                    if (s == 0)
                        record.strand = '+';
                    else
                        record.strand = '-';

                    unsigned prev_cs = t;
                    double scoresSum = (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) );
                    std::stringstream ss_indivScores;
                    ss_indivScores << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) ) << ';';
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
                            double secondBest = 0.0;
                            for (unsigned k = 0; k < 4; ++k)
                            {
                                if (k != (unsigned)data.states[s][i][t] && data.statePosteriors[s][k][i][t] > secondBest)
                                    secondBest = data.statePosteriors[s][k][i][t];
                            }                   

                            scoresSum += (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) );
                            ss_indivScores << (double)log(data.statePosteriors[s][data.states[s][i][t]][i][t] / std::max(secondBest, DBL_MIN) ) << ';';
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
                    writeRecord(outBed, record);
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

#endif
