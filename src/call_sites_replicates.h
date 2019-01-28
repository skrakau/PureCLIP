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



#ifndef APPS_HMMS_CALL_SITES_REPLICATES_H_
#define APPS_HMMS_CALL_SITES_REPLICATES_H_

#define HMM_PARALLEL 1

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <unistd.h>     
#include <sys/stat.h>
#include <errno.h>
#include <seqan/modifier.h>
#include <seqan/bed_io.h>

#include "parse_alignments.h"
#include "hmm_1.h"
#include "prepro_mle.h"

#include "density_functions_reg.h"
#include "density_functions_crosslink.h"
#include "density_functions_crosslink_reg.h"

using namespace seqan;

//         Infix<String<__uint16> >::Type truncCounts;
//         unsigned contigId; 
// 
//         String<__uint32>    nEstimates;      
//         String<double>      kdes;       // used for 'enriched'. 'non-enriched' classification
//         String<double>      kdesN;    // used to estimate the binomial n parameters (decoupled, might be useful e.g. for longer crosslink clusters) 
//         String<double>      rpkms;      // TODO change name -> e.g. bgSignal
//         String<float>       fimoScores; // for each t: one motif score
//         String<char>        motifIds; // for each t: one motif score



void copyAndClipData(Data &newData, Data &oldData, unsigned s, unsigned i, unsigned beginPos, unsigned endPos, unsigned newSetPos)
{
//     std::cout << "  i " << i << std::endl;
//     std::cout << "  beginPos " << beginPos << std::endl;
//     std::cout << "  endPos " << endPos << std::endl;
//     std::cout << "  newSetPos " << newSetPos << std::endl;
    
    Observations obs;          
    obs.truncCounts = infix(oldData.setObs[s][i].truncCounts, beginPos, endPos);    // infix of infix

//     std::cerr << "length(oldData.setObs[" << s << "][ " << i << "].truncCounts)=" << length(oldData.setObs[s][i].truncCounts) << std::endl;
//     std::cerr << "length(host(oldData.setObs[s][i].truncCounts))=" << length(host(oldData.setObs[s][i].truncCounts)) << std::endl;
//     std::cerr << "beginPosition(oldData.setObs[s][i].truncCounts)=" << beginPosition(oldData.setObs[s][i].truncCounts) << std::endl;
//     std::cerr << "endPosition(oldData.setObs[s][i].truncCounts)=" << endPosition(oldData.setObs[s][i].truncCounts) << std::endl;
//     std::cerr << std::endl;
//     std::cerr << "length(obs.truncCounts)=" << length(obs.truncCounts) << std::endl;
//     std::cerr << "length(host(obs.truncCounts))=" << length(host(obs.truncCounts)) << std::endl;
//     std::cerr << "beginPosition(obs.truncCounts)=" << beginPosition(obs.truncCounts) << std::endl;
//     std::cerr << "endPosition(obs.truncCounts)=" << endPosition(obs.truncCounts) << std::endl;
//     std::cerr << std::endl;
//     std::cerr << "host(oldData.setObs[s][i].truncCounts)[beginPosition(oldData.setObs[s][i].truncCounts)]=" << host(oldData.setObs[s][i].truncCounts)[beginPosition(oldData.setObs[s][i].truncCounts)] << std::endl;
// 
//     std::cerr << "host(obs.truncCounts)[beginPosition(obs.truncCounts)]=" << host(obs.truncCounts)[beginPosition(obs.truncCounts)] << std::endl;
//     std::cerr << "infix(host(obs.truncCounts), beginPosition(obs.truncCounts), endPosition(obs.truncCounts))[0]=" << infix(host(obs.truncCounts), beginPosition(obs.truncCounts), endPosition(obs.truncCounts))[0] << std::endl;
//     std::cerr << std::endl;
//     std::cerr << "obs.truncCounts[0]=" << obs.truncCounts[0] << std::endl;

    obs.nEstimates = infix(oldData.setObs[s][i].nEstimates, beginPos, endPos);
    obs.kdes = infix(oldData.setObs[s][i].kdes, beginPos, endPos);
    obs.kdesN = infix(oldData.setObs[s][i].kdesN, beginPos, endPos);
    if (!empty(oldData.setObs[s][i].rpkms)) obs.rpkms = infix(oldData.setObs[s][i].rpkms, beginPos, endPos);
    if (!empty(oldData.setObs[s][i].fimoScores)) obs.fimoScores = infix(oldData.setObs[s][i].fimoScores, beginPos, endPos);
     if (!empty(oldData.setObs[s][i].motifIds)) obs.motifIds = infix(oldData.setObs[s][i].motifIds, beginPos, endPos);
    obs.discard = false;
    obs.contigId = oldData.setObs[s][i].contigId;

    appendValue(newData.setObs[s], obs);  
    appendValue(newData.setPos[s], newSetPos);  

    // state posteriors: still empty
    // states: still empty
}





// TODO TEST !!!
bool intersect_replicateIntervals(String<Data> &data_replicates)
{
    // currently only for 2 replicates!?
    if (length(data_replicates) > 2) 
    {
        std::cerr << "ERROR: not supporting > 2 replicates yet!" << std::endl;
        return false;
    }

    Data newData_rep1;
    resize(newData_rep1.setObs, 2);
    resize(newData_rep1.setPos, 2);
    resize(newData_rep1.statePosteriors, 2);
    resize(newData_rep1.states, 2); 
    Data newData_rep2 = newData_rep1;

    // intersect & clip
    for (unsigned s = 0; s < 2; ++s)
    {
        // FORWARD (&& REVERSE ?)
        unsigned it1 = 0;
        unsigned it2 = 0;
        while (it1 < length(data_replicates[0].setObs[s]) && it2 < length(data_replicates[1].setObs[s]))
        {
            // discard not intersecting or discarded!!!
            if ((data_replicates[0].setPos[s][it1] + data_replicates[0].setObs[s][it1].length() <= data_replicates[1].setPos[s][it2]) ||
                 data_replicates[0].setObs[s][it1].discard) 
            {
//                 if (s == 0 && data_replicates[0].setPos[s][it1]  <= 46912654 && data_replicates[0].setPos[s][it1] >= 46800000)
//                     std::cout << "Found no intersection for rep1 interval .." <<  data_replicates[0].setPos[s][it1] << " - " << (data_replicates[0].setPos[s][it1] + data_replicates[0].setObs[s][it1].length()) << std::endl;
                
                ++it1;
            }
            else if ((data_replicates[1].setPos[s][it2] + data_replicates[1].setObs[s][it2].length() <= data_replicates[0].setPos[s][it1]) ||
                      data_replicates[1].setObs[s][it2].discard)
            {
//                 if (s == 0 && data_replicates[1].setPos[s][it2]  <= 46912654 && data_replicates[1].setPos[s][it2] >= 46800000)
//                     std::cout << "Found no intersection for rep2 interval .." <<  data_replicates[1].setPos[s][it2] << " - " << (data_replicates[1].setPos[s][it2] + data_replicates[1].setObs[s][it2].length()) << std::endl;

                ++it2;
            }
            else
            {
//                 if (s == 0 && data_replicates[0].setPos[s][it1]  <= 46912654 && data_replicates[0].setPos[s][it1] >= 46800000)
//                 {
//                     std::cout << "Found intersection .." << it1 << " " << it2 << std::endl;
//                     std::cout << "  rep1: " << data_replicates[0].setPos[s][it1] << " - " << (data_replicates[0].setPos[s][it1] + data_replicates[0].setObs[s][it1].length()) << std::endl;
//                     std::cout << "  rep2: " << data_replicates[1].setPos[s][it2] << " - " << (data_replicates[1].setPos[s][it2] + data_replicates[1].setObs[s][it2].length()) << std::endl;
//                 }
                
                unsigned rep1_frontClipPos = 0;
                unsigned rep2_frontClipPos = 0;
                unsigned rep1_backClipPos = data_replicates[0].setObs[s][it1].length();
                unsigned rep2_backClipPos = data_replicates[1].setObs[s][it2].length();
                unsigned new_setPos = 0;

                if (data_replicates[0].setPos[s][it1] < data_replicates[1].setPos[s][it2])  // clip rep1 interval front 
                {
                    rep1_frontClipPos = data_replicates[1].setPos[s][it2] - data_replicates[0].setPos[s][it1];  
                    new_setPos = data_replicates[1].setPos[s][it2];
                }
                else                                                                     // clip rep2 interval front
                {
                    rep2_frontClipPos = data_replicates[0].setPos[s][it1] - data_replicates[1].setPos[s][it2];
                    new_setPos = data_replicates[0].setPos[s][it1]; 
                }

                if (data_replicates[0].setPos[s][it1] + data_replicates[0].setObs[s][it1].length() < data_replicates[1].setPos[s][it2] + data_replicates[1].setObs[s][it2].length())    // clip rep2 interval back
                {
                    rep2_backClipPos = data_replicates[0].setPos[s][it1] + data_replicates[0].setObs[s][it1].length() - data_replicates[1].setPos[s][it2];
                }
                else                                                                                                                                                                // clip rep1 interval back
                {
                    rep1_backClipPos = data_replicates[1].setPos[s][it2] + data_replicates[1].setObs[s][it2].length() - data_replicates[0].setPos[s][it1];
                }
                
//                 if (s == 0 && data_replicates[0].setPos[s][it1]  <= 46912654 && data_replicates[0].setPos[s][it1] >= 46800000)
//                 {
//                     std::cout << "Copy and clip .." << std::endl;
//                     std::cout << "  rep1: new start: " << new_setPos << " clip at: " << rep1_frontClipPos << " - " << rep1_backClipPos << std::endl;
//                     std::cout << "  rep2: new start: " << new_setPos << " clip at: " << rep2_frontClipPos << " - " << rep2_backClipPos << std::endl;
//                 }
                
                copyAndClipData(newData_rep1, data_replicates[0], s, it1, rep1_frontClipPos, rep1_backClipPos, new_setPos);
                copyAndClipData(newData_rep2, data_replicates[1], s, it2, rep2_frontClipPos, rep2_backClipPos, new_setPos);

                if (data_replicates[0].setPos[s][it1] + data_replicates[0].setObs[s][it1].length() < data_replicates[1].setPos[s][it2] + data_replicates[1].setObs[s][it2].length()) 
                    ++it1;
                else
                    ++it2;
            }
        }
    }
    // check rep 0 & it1, rep 1 & it2
    
    data_replicates[0] = newData_rep1;
    data_replicates[1] = newData_rep2;

    return true;
}



template<typename TDOUBLE>
long double getPseudo_eProb(Observations &obs, unsigned t, ZTBIN<TDOUBLE> &bin2, AppOptions &options)
{
    return bin2.getDensity(1, 2*obs.nEstimates[t], options); 
}

template<typename TDOUBLE>
long double getPseudo_eProb(Observations &obs, unsigned t, ZTBIN_REG<TDOUBLE> &bin2, AppOptions &options)
{
    unsigned mId = obs.motifIds[t];
    long double bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*obs.fimoScores[t]));

    return bin2.getDensity(1, 2*obs.nEstimates[t], bin2_pred, options);
}



template<typename TGAMMA, typename TBIN, typename TDOUBLE>
HMM<TGAMMA, TBIN, TDOUBLE> merge_HMMs(String<HMM<TGAMMA, TBIN, TDOUBLE> > &hmms_replicates, String<TBIN> &bin2, AppOptions &options)
{
    HMM<TGAMMA, TBIN, TDOUBLE> mergedHmm = hmms_replicates[0];
//    HMM<TGAMMA, TBIN, TDOUBLE>mergedHmm(4, data_replicates[rep].setObs, data_replicates[rep].setPos, contigLen); //hmms_replicates[0];
 
    std::cout << "   set initial probs" << std::endl;        
    // init probs
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(mergedHmm.setObs[s]); ++i)
        {
            for (unsigned k = 0; k < 4; ++k)
            {
                double sum = 0.0;
                for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                    sum += hmms_replicates[rep].initProbs[s][i][k];                     // TODO use position-specific post. probs. instead?

                mergedHmm.initProbs[s][i][k] = sum/(double)length(hmms_replicates);
            }
        }
    }
    std::cout << "   set trans probs" << std::endl;            
    // trans matrix
    for (unsigned k_1 = 0; k_1 < 4; ++k_1)
    {
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
        {
            double sum = 0.0;
            for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                sum += hmms_replicates[rep].transMatrix[k_1][k_2];

            mergedHmm.transMatrix[k_1][k_2] = sum/(double)(double)length(hmms_replicates);
        }
    }
    std::cout << "   set emission probs" << std::endl;            
    /////////////////////
    // compute new eProbs
    /////////////////////
    // pseudo-counbt -> pseudo-eprob ?
    // avoid becoming zero, if no read start in one replicate
    // n -> k=1, 2*n
    // we do not have 'crosslinkin' eprobs here, only joint eprobs!
    // if k = 0,
    // TODO pseudo-eprobs also for non-crosslink state!
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(mergedHmm.setObs[s]); ++i)
        {
            for (unsigned t = 0; t < mergedHmm.setObs[s][i].length(); ++t)        // length(data.states[s][i])
            {
                for (unsigned k = 0; k < 3; ++k)
                {
                    mergedHmm.eProbs[s][i][t][k] = 1.0; 
                    for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                        mergedHmm.eProbs[s][i][t][k] *= hmms_replicates[rep].eProbs[s][i][t][k];
                }
                // for crosslink state
                mergedHmm.eProbs[s][i][t][3] = 1.0; //hmms_replicates[0].eProbs[s][i][t][k];
                for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                {
                    if (options.use_pseudoEProb && hmms_replicates[rep].setObs[s][i].truncCounts[t] == 0)  // add pseudo-count
                        mergedHmm.eProbs[s][i][t][3] *= hmms_replicates[rep].eProbs[s][i][t][2] * getPseudo_eProb(hmms_replicates[rep].setObs[s][i], t, bin2[rep], options);
                    else
                        mergedHmm.eProbs[s][i][t][3] *= hmms_replicates[rep].eProbs[s][i][t][3];
                }
            }
        }
    }
    // state post. probs. -> will be overwritten anyways
    
    std::cout << "   set setObs" << std::endl;        
    // estimated n (used for threshold when learning trans. probs 2-> 3)
    // other obs ?
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(mergedHmm.setObs[s]); ++i)
        {
            for (unsigned t = 0; t < mergedHmm.setObs[s][i].length(); ++t)        // length(data.states[s][i])
            {
                unsigned minN = 10000;
                for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                    if (hmms_replicates[rep].setObs[s][i].nEstimates[t] < minN) minN = hmms_replicates[rep].setObs[s][i].nEstimates[t];
                
                mergedHmm.setObs[s][i].nEstimates[t] = minN;

                // set others to 0 for now
                //mergedHmm.setObs[s][i].truncCounts[t] = 0;        // reference, cannot be changed!        
                mergedHmm.setObs[s][i].kdes[t] = 0.0;
                mergedHmm.setObs[s][i].kdesN[t] = 0.0;    
            }
        }
    }
    return mergedHmm;
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE> 
bool HMM<TGAMMA, TBIN, TDOUBLE>::updateTransProbsAndPostDecode(AppOptions &options)
{
    double old_transProb02 = this->transMatrix[0][2];
    double old_transProb23 = this->transMatrix[2][3];
    for (unsigned iter = 0; iter < options.maxIter_bw; ++iter)
    {
        std::cout << ".. " << iter << "th iteration " << std::endl;
  
        std::cout << "                        computeStatePosteriorsFB() " << std::endl;
        if (!computeStatePosteriorsFBupdateTrans(options))      // also updates initial probabilities
            return false;
       
        std::cout << "*** Transition probabilitites ***" << std::endl;
        for (unsigned k_1 = 0; k_1 < this->K; ++k_1)
        {
            std::cout << "    " << k_1 << ": ";
            for (unsigned k_2 = 0; k_2 < this->K; ++k_2)
                std::cout << this->transMatrix[k_1][k_2] << "  ";
            std::cout << std::endl;
        }
        std::cout << std::endl;

        if (std::fabs(this->transMatrix[2][3] - old_transProb23) > 0.00001 && std::fabs(this->transMatrix[0][2] - old_transProb02) > 0.00001)
        {
            std::cout << " **** Convergence ! **** " << std::endl;
            break;
        }
        old_transProb02 = this->transMatrix[0][2];
        old_transProb23 = this->transMatrix[2][3];
    }
    return true;
}



// in parallel for replicates
template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool applyHMM(Data &newData,
              String<Data> &data_replicates, 
              String<String<String<long double> > > &transMatrices_1,
              String<TGAMMA> &d1, String<TGAMMA> &d2, String<TBIN> &bin1, String<TBIN> &bin2,
              TDOUBLE /**/,
              unsigned &contigLen,
              AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    String<HMM<TGAMMA, TBIN, TDOUBLE> > hmms_replicates;
    reserve(hmms_replicates, length(data_replicates));
    for (unsigned rep = 0; rep < length(data_replicates); ++rep)
    {
        if (options.verbosity >= 1) std::cout << "Apply learned HMM for replicate: " << rep << std::endl;
        appendValue(hmms_replicates, HMM<TGAMMA, TBIN, TDOUBLE>(4, data_replicates[rep].setObs, data_replicates[rep].setPos, contigLen));
        auto & hmm = hmms_replicates[rep];
        hmm.transMatrix = transMatrices_1[rep];
        std::cout << "   applyParameters" << std::endl;        
        if (!hmm.applyParameters(d1[rep], d2[rep], bin1[rep], bin2[rep], options))
             return false;

//         HMM<TGAMMA, TBIN, TDOUBLE>hmm(4, data_replicates[rep].setObs, data_replicates[rep].setPos, contigLen);
//         std::cout << "   set transMatrix" << std::endl;
//         hmm.transMatrix = transMatrices_1[rep];
//         std::cout << "   applyParameters" << std::endl;
//         if (!hmm.applyParameters(d1[rep], d2[rep], bin1[rep], bin2[rep], options))            
//             return false;
    }

    // merge HMMs and multiply eProbs!
    if (options.verbosity >= 1) std::cout << "Merge HMMs from replicates ... " << std::endl;
    HMM<TGAMMA, TBIN, TDOUBLE> mergedHmm = merge_HMMs(hmms_replicates, bin2, options);

    // update
    if (options.verbosity >= 1) std::cout << "Update transition probabilities and get final posterior probabilities ... " << std::endl;    
    mergedHmm.updateTransProbsAndPostDecode(options);

    // get states
    if (options.posteriorDecoding)
        mergedHmm.posteriorDecoding(newData.states);
    else
        mergedHmm.viterbi(newData.states);

//     if (options.useCov_RPKM)    // NOTE: otherwise not necessary, since gamma1.k <= 1
//         mergedHmm.rmBoarderArtifacts(newData.states, d1[0]);    // TODO 

    newData.statePosteriors = mergedHmm.statePosteriors;
    newData.setObs = mergedHmm.setObs;
    newData.setPos = mergedHmm.setPos;

#ifdef HMM_PROFILE
    Times::instance().time_applyHMM += (sysTime() - timeStamp);
#endif

    return true;
}


template <typename TGamma, typename TBIN, typename TDOUBLE, typename TOptions>
bool doIt(String<TGamma> &gamma1, String<TGamma> &gamma2, String<TBIN> &bin1, String<TBIN> &bin2, TDOUBLE /**/, TOptions &options)
{
#if HMM_PARALLEL
    omp_set_num_threads(options.numThreads);
#endif  

#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    typedef  FragmentStore<>    TStore;
    TStore store;
    if (options.verbosity >= 1) std::cout << "Loading reference ... " << std::endl;
    
    try {
        if (!loadContigs(store, toCString(options.refFileName)))
        {
            std::cerr << "ERROR: Can't load reference sequence from file '" << options.refFileName << "'" << std::endl;
            return 1;
        }
    } catch (std::exception &e){
        std::cerr << "ERROR: Can't load reference sequence from file '" << options.refFileName << "': " << e.what() << std::endl;
        return 1;
    }
    loadIntervals(options, store);
    loadApplyChrs(options, store);
    if (options.numThreadsA == 0)
    {
        options.numThreadsA = std::min((int)options.numThreads, (int)length(options.intervals_contigIds));
        if (options.verbosity >= 2)  std::cout << "Set no. of threads for applying HMM to: " << options.numThreadsA << std::endl;
    }
 
#if SEQAN_HAS_ZLIB
    if (options.verbosity > 1) std::cout << "SEQAN_HAS_ZLIB" << std::endl;
#else
    std::cout << "WARNING: zlib not available !" << std::endl;
#endif

    // ******************  set some parameters
    // some precision related:
    options.min_nligf = std::nextafter((long double)1.0, (long double)0.0);
    if (options.verbosity >= 2) std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << "Max. nligf: " << options.min_nligf << std::setprecision(6) << std::endl;
    int db_min_exp = DBL_MIN_10_EXP;
    int ldb_min_exp = LDBL_MIN_10_EXP;
    std::cout << "DBL_MIN_10_EXP: " << db_min_exp << " LDBL_MIN_10_EXP: " << ldb_min_exp << std::endl;
 
    // some other thresholds and setings:
    if (options.binSize == 0.0) options.binSize = options.bandwidth * 2; 
    options.intervalOffset = options.bandwidth * 2;  
    options.prior_kdeThreshold = options.prior_enrichmentThreshold * getGaussianKernelDensity(0.0/(double)options.bandwidth)/(double)options.bandwidth;
    if (options.verbosity >= 1) std::cout << " computed prior_kdeThreshold: " << (options.prior_enrichmentThreshold * getGaussianKernelDensity(0.0/(double)options.bandwidth)/(double)options.bandwidth) << std::endl;

    if (options.useKdeThreshold == 0.0 && !options.useCov_RPKM) // TODO use boolean user option
        options.useKdeThreshold = getGaussianKernelDensity(0.0/(double)options.bandwidth)/(double)options.bandwidth + 0.0001;  // corresponds to KDE value at singleton read start 
    else if (options.useKdeThreshold == 0.0 && options.useCov_RPKM)
    {
        options.useKdeThreshold = getGaussianKernelDensity(0.0/(double)options.bandwidth)/(double)options.bandwidth; 
        if (options.mrtf_kdeSglt)
            options.minRPKMtoFit = log(options.useKdeThreshold) + 0.0001;

    }
    if (options.verbosity >= 1) std::cout << "Use bandwidth: " << options.bandwidth << std::endl;
    if (options.verbosity >= 1) std::cout << "Use KDE threshold: " << options.useKdeThreshold << std::endl;
    if (options.verbosity >= 1) std::cout << "Use bandwidth to estimate n: " << options.bandwidthN << std::endl;

    // if required, determine n threshold for learning of binomial parameters and trans. probs (2-> 2/3)
    if (options.get_nThreshold)
    {
        // require at least a mean of 2 read start counts for 'crosslink' state
        options.nThresholdForTransP = ceil(2.0/options.p2);
        options.nThresholdForP = options.nThresholdForTransP;
        std::cout << "Set n threshold used for learning of binomial p parameters and transition probabilities '2' -> '2'/'3' to: " << options.nThresholdForTransP << std::endl;
    } 


    // Read control BAI index.
    BamIndex<Bai> inputBaiIndex;
    if (options.useCov_RPKM && !empty(options.inputBaiFileName))
    {
        if (!open(inputBaiIndex, toCString(options.inputBaiFileName)))
        {
            std::cerr << "ERROR: Could not read input control BAI index file " << options.inputBaiFileName << "\n";
            return 1;
        }
    }

    ////////////////////////////////////////////////////
    // for each replicate, learn HMM parameters
    ////////////////////////////////////////////////////
    // gamma1, gamma2, bin1, bin2 
    // TODO should be given as list already
    String<double> slr_NfromKDE_b0;
    String<double> slr_NfromKDE_b1;
    resize(slr_NfromKDE_b0, length(options.baiFileNames), 0.0);
    resize(slr_NfromKDE_b1, length(options.baiFileNames), 0.0);
    String<String<String<long double> > > transMatrices;
    resize(transMatrices, length(options.baiFileNames));
    String<BamIndex<Bai> > baiIndices;
    resize(baiIndices, length(options.baiFileNames));

    for (unsigned rep = 0; rep < length(options.baiFileNames); ++rep)
    {
        if (options.verbosity >= 1) std::cout << "Learn HMM parameters for replicate: " << rep << std::endl;
        
        // Read BAI index.
        BamIndex<Bai> baiIndex;
        if (!open(baiIndex, toCString(options.baiFileNames[rep])))
        {
            std::cerr << "ERROR: Could not read BAI index file " << options.baiFileNames[rep] << "\n";
            return 1;
        }
        baiIndices[rep] = baiIndex;

        // *****************

        String<ContigObservations> contigObservationsF;
        String<ContigObservations> contigObservationsR;
        resize(contigObservationsF, length(options.intervals_contigIds), Exact());
        resize(contigObservationsR, length(options.intervals_contigIds), Exact());

        Data data;
        resize(data.setObs, 2);
        resize(data.setPos, 2);
        resize(data.statePosteriors, 2);
        resize(data.states, 2);
        bool stop = false;
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
            for (unsigned i = 0; i < length(options.intervals_contigIds); ++i)
            {
                unsigned contigId = options.intervals_contigIds[i];

                int r = loadObservations(contigObservationsF[i], contigObservationsR[i], contigId, options.bamFileNames[rep], baiIndices[rep], store, options);
                if (r == 1)
                {
                    stop = true; 
                }
                else if (r == 0)
                {
                    String<double> contigCovsF;
                    String<double> contigCovsR;
                    loadCovariates(contigCovsF, contigCovsR, contigId, store, options); 
                    String<String<float> > contigCovsFimo;
                    String<String<char> > motifIds;
                    loadMotifCovariates(contigCovsFimo, motifIds, contigId, store, options); 

                    // Extract covered intervals for learning
                    unsigned i1 = options.intervals_positions[i][0];    // interval begin
                    unsigned i2 = options.intervals_positions[i][1];    // interval end
                    Data c_data;                
                    resize(c_data.setObs, 2);
                    resize(c_data.setPos, 2);
                    resize(c_data.statePosteriors, 2);
                    resize(c_data.states, 2);

                    extractCoveredIntervals(c_data, contigObservationsF[i], contigObservationsR[i], contigCovsF, contigCovsR, contigCovsFimo, motifIds, contigId, i1, i2, options.excludePolyAFromLearning, options.excludePolyTFromLearning, true, store, options); 

                    SEQAN_OMP_PRAGMA(critical)
                        append(data, c_data);  
                }
            }
        if (stop) return 1;

        // learn KDE - N relationship on all contigs used for other parameter learning as well
        for (unsigned s = 0; s < 2; ++s)
            for (unsigned i = 0; i < length(data.setObs[s]); ++i)
                data.setObs[s][i].computeKDEs(options);

        if (options.estimateNfromKdes) 
            computeSLR(slr_NfromKDE_b0[rep], slr_NfromKDE_b1[rep], data, options);

        // precompute KDE values, estimate Ns, etc.
        preproCoveredIntervals(data, slr_NfromKDE_b0[rep], slr_NfromKDE_b1[rep], inputBaiIndex, store, true, options);

        gamma1[rep].tp = options.useKdeThreshold;
        gamma2[rep].tp = options.useKdeThreshold;       // left tuncated    

        if (options.verbosity >= 1) std::cout << "Prior ML estimation of density distribution parameters using predefined cutoff ..." << std::endl;
        prior_mle(gamma1[rep], gamma2[rep], data, options);
        estimateTransitions(transMatrices[rep], gamma1[rep], gamma2[rep], bin1[rep], bin2[rep], data, options);

        unsigned contigLen = 0; // should not be used within learning
        if (!learnHMM(data, transMatrices[rep], gamma1[rep], gamma2[rep], bin1[rep], bin2[rep], (TDOUBLE)0.0, contigLen, options))
            return 1;

        clear(contigObservationsF);
        clear(contigObservationsR);
        clear(data);
    }


    if (options.verbosity >= 1) std::cout << "Apply learned parameters to whole genome  ..." << std::endl;
#if HMM_PARALLEL
        omp_set_num_threads(options.numThreadsA);
#endif  


    // store contig-wise bedRecords
    String<String<BedRecord<Bed6> > > bedRecords_sites;
    String<String<BedRecord<Bed6> > > bedRecords_regions;
    resize(bedRecords_sites, length(options.applyChr_contigIds), Exact());
    resize(bedRecords_regions, length(options.applyChr_contigIds), Exact());

#ifdef HMM_PROFILE
    double timeStamp2 = sysTime();
#endif

    for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
    {
        unsigned contigId = options.applyChr_contigIds[i];
        unsigned contigLen = length(store.contigStore[contigId].seq);

        if (options.verbosity >= 1) std::cout << "  " << store.contigNameStore[contigId] << std::endl;

        String<ContigObservations> contigObservationsF;
        String<ContigObservations> contigObservationsR;
        resize(contigObservationsF, length(options.baiFileNames));
        resize(contigObservationsR, length(options.baiFileNames));
        String<Data> data_replicates;
        resize(data_replicates, length(options.baiFileNames));    

        bool stop = false;
//#if HMM_PARALLEL
//    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreadsA)) 
//#endif 
        for (unsigned rep = 0; rep < length(options.baiFileNames); ++rep)
        {
            if (options.verbosity >= 1) std::cout << "Get preprocessed intervals for replicate: " << rep << std::endl;
            
            // for each replicate get preprocessed covered intervals
            int r = loadObservations(contigObservationsF[rep], contigObservationsR[rep], contigId, options.bamFileNames[rep],baiIndices[rep], store, options);
            if (r == 1)
            {
                return 1;
                //stop = true; 
            }
            else if (r == 0)
            {
                String<double> c_contigCovsF;
                String<double> c_contigCovsR;
                loadCovariates(c_contigCovsF, c_contigCovsR, contigId, store, options); 
                String<String<float> > c_contigCovsFimo;
                String<String<char> > c_motifIds;
                loadMotifCovariates(c_contigCovsFimo, c_motifIds, contigId, store, options); 

                // Extract covered intervals
                unsigned i1 = 0;    
                unsigned i2 = length(store.contigStore[contigId].seq);    
                Data c_data;                
                resize(c_data.setObs, 2);
                resize(c_data.setPos, 2);
                resize(c_data.statePosteriors, 2);
                resize(c_data.states, 2); 
                extractCoveredIntervals(c_data, contigObservationsF[rep], contigObservationsR[rep], c_contigCovsF, c_contigCovsR, c_contigCovsFimo, c_motifIds, contigId, i1, i2, options.excludePolyA, options.excludePolyT, false, store, options); 

                if (!empty(c_data.setObs[0]) || !empty(c_data.setObs[1]))   
                {
                    preproCoveredIntervals(c_data, slr_NfromKDE_b0[rep], slr_NfromKDE_b1[rep], inputBaiIndex, store, false, options);
                    data_replicates[rep] = c_data;
                }
            }
        }
        //if (stop) return 1;

        // get intersection of covered intervals and clip individual observations
        if (options.verbosity >= 1) std::cout << "Intersect covered intervals from replicates ... " << std::endl;
        intersect_replicateIntervals(data_replicates);
        
        Data repData = data_replicates[0];
        // build individual HMMs on clipped intervals, merge HMMs, update trans. probs, compute post. probs and get states
        if (!applyHMM(repData, data_replicates, transMatrices, gamma1, gamma2, bin1, bin2, (TDOUBLE)0.0, contigLen, options))
        {
            return 1;
        }

        writeStates(bedRecords_sites[i], repData, store, contigId, options);  
        if (!empty(options.outRegionsFileName))
        {
            writeRegions(bedRecords_regions[i], repData, store, contigId, options);              
        }
    }

    if (options.verbosity >= 2) std::cout << "Write bedRecords to BED file ... " << options.outFileName << std::endl;
    // crosslink sites
    BedFileOut outBed(toCString(options.outFileName)); 
    for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
    {
        for (unsigned j = 0; j < length(bedRecords_sites[i]); ++j)
            writeRecord(outBed, bedRecords_sites[i][j]);
    }
    // binding regions
    if (!empty(options.outRegionsFileName))
    {
        BedFileOut outBed2(toCString(options.outRegionsFileName)); 
        for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
        {
            for (unsigned j = 0; j < length(bedRecords_regions[i]); ++j)
                writeRecord(outBed2, bedRecords_regions[i][j]);
        }
    }

#ifdef HMM_PROFILE
    Times::instance().time_applyHMM2 = sysTime() - timeStamp2;
#endif

#ifdef HMM_PROFILE
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
    std::cout << "  Time needed for loadObservations: " << Times::instance().time_loadObservations/60.0 << "min" << std::endl;
    std::cout << "  Time needed for learnHMM: " << Times::instance().time_learnHMM/60.0 << "min" << std::endl;
    std::cout << "  Time needed for applyHMM: " << Times::instance().time_applyHMM/60.0 << "min" << std::endl;
    //std::cout << "  Time needed for applyHMM2: " << Times::instance().time_applyHMM2/60.0 << "min" << std::endl;
#endif


    CharString fileNameParams;
    if (!empty(options.parFileName))
        fileNameParams = options.parFileName;
    else
    {
        fileNameParams = prefix(options.outFileName, length(options.outFileName)-4);
        append(fileNameParams, ".params");
    }
    std::ofstream out(toCString(fileNameParams), std::ios::binary | std::ios::out);
    out << "options.useKdeThreshold" << '\t' << options.useKdeThreshold << std::endl;
    out << std::endl;
    for (unsigned rep = 0; rep < length(options.baiFileNames); ++rep)
    {
        out << "Parameter learned for replicate: " << rep << std::endl;
        printParams(out, gamma1[rep], 1);
        printParams(out, gamma2[rep], 2);
        printParams(out, bin1[rep], 1);
        printParams(out, bin2[rep], 2);
        printParams(out, transMatrices[rep]);
        out << std::endl;
        out << "slr_NfromKDE.b0" << '\t' << slr_NfromKDE_b0[rep] << std::endl;
        out << "slr_NfromKDE.b1" << '\t' << slr_NfromKDE_b1[rep] << std::endl;  
    }

    return 0;
}

#endif
