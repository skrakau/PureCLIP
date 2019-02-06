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
            // discard: not intersecting or discarded!!!
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



// long double getPseudo_eProb(Observations &obs, unsigned t, ZTBIN &bin2, AppOptions &options)
// {
//     return bin2.getDensity(1, 2*obs.nEstimates[t], options); 
// }
// 
// long double getPseudo_eProb(Observations &obs, unsigned t, ZTBIN_REG &bin2, AppOptions &options)
// {
//     unsigned mId = obs.motifIds[t];
//     long double bin2_pred = 1.0/(1.0+exp(-bin2.b0 - bin2.regCoeffs[mId]*obs.fimoScores[t]));
// 
//     return bin2.getDensity(1, 2*obs.nEstimates[t], bin2_pred, options);
// }



template<typename TGAMMA, typename TBIN>
HMM<TGAMMA, TBIN> merge_HMMs(String<HMM<TGAMMA, TBIN> > &hmms_replicates, String<ModelParams<TGAMMA, TBIN> > &modelParams, AppOptions &options)
{
    HMM<TGAMMA, TBIN> mergedHmm = hmms_replicates[0];

    if (length(options.baiFileNames) > 1)
    {
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
        // pseudo-count -> pseudo-eprob ?
        // avoid becoming zero, if no read start in one replicate
        // n -> k=1, 2*n
        // we do not have 'crosslinkin' eprobs here, only joint eprobs!
        // if k = 0,
        for (unsigned s = 0; s < 2; ++s)
        {
            for (unsigned i = 0; i < length(mergedHmm.setObs[s]); ++i)
            {
                for (unsigned t = 0; t < mergedHmm.setObs[s][i].length(); ++t)
                {
                    for (unsigned k = 0; k < 3; ++k)
                    {
                        mergedHmm.eProbs[s][i][t][k] = 0.0; 
                        for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                            mergedHmm.eProbs[s][i][t][k] += hmms_replicates[rep].eProbs[s][i][t][k];    // log-space
                    }
                    // for crosslink state
                    mergedHmm.eProbs[s][i][t][3] = 0.0;
                    for (unsigned rep = 0; rep < length(hmms_replicates); ++rep)
                    {
//                         if (options.use_pseudoEProb && hmms_replicates[rep].setObs[s][i].truncCounts[t] == 0)  // add pseudo-count
//                             mergedHmm.eProbs[s][i][t][3] += hmms_replicates[rep].eProbs[s][i][t][2] * getPseudo_eProb(hmms_replicates[rep].setObs[s][i], t, modelParams[rep].bin2, options);
//                         else
                        mergedHmm.eProbs[s][i][t][3] += hmms_replicates[rep].eProbs[s][i][t][3];    // log-space
                    }
                }
            }
        }
        // state post. probs. -> will be overwritten anyways

        if (options.verbosity >= 1) std::cout << "   set setObs" << std::endl;        
        // estimated n (used for threshold when learning trans. probs 2-> 3)
        // other obs not used anymore
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
    }
    return mergedHmm;
}


template<typename TGAMMA, typename TBIN> 
bool HMM<TGAMMA, TBIN>::updateTransAndPostProbs(AppOptions &options)
{
    double old_transProb02 = this->transMatrix[0][2];
    double old_transProb23 = this->transMatrix[2][3];
    for (unsigned iter = 0; iter < options.maxIter_bw; ++iter)
    {
        std::cout << ".. " << iter << "th iteration " << std::endl;
  
        std::cout << "                        computeStatePosteriorsFB() " << std::endl;
        if (!computeStatePosteriorsFBupdateTrans(options))      // also updates initial probabilities
            return false;
       
        std::cout << "*** Transition probabilities ***" << std::endl;
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



/////////////////////////////////////////////
// TEST
void setUp(String<Data> &data_replicates)
{
    resize(data_replicates, 2);

    Data newData_rep1;
    resize(newData_rep1.setObs, 2);
    resize(newData_rep1.setPos, 2);
    resize(newData_rep1.statePosteriors, 2);
    resize(newData_rep1.states, 2);
    Data newData_rep2 = newData_rep1;

    unsigned no_t = 10;
    Observations obsShort;
    String<__uint16> truncs;
    resize(truncs, 100, 0);    
    obsShort.truncCounts = infix(truncs, 10, 20);
    resize(obsShort.nEstimates, no_t, 0);
    resize(obsShort.kdes, no_t, 0.0);
    resize(obsShort.kdesN, no_t, 0.0);
    obsShort.discard = false;
    obsShort.contigId = 0;

    no_t = 20;
    Observations obsLong;
    obsLong.truncCounts = infix(truncs, 10, 30);
    resize(obsLong.nEstimates, no_t, 0);
    resize(obsLong.kdes, no_t, 0.0);
    resize(obsLong.kdesN, no_t, 0.0);
    obsLong.discard = false;
    obsLong.contigId = 0;

    // ---    ---   ---        ----     ---  
    //  ---        -----  ---   --    ---

    for (unsigned s = 0; s < 2; ++s)
    {
        resize(newData_rep1.setObs[s], 5);
        resize(newData_rep1.setPos[s], 5);  
        resize(newData_rep2.setObs[s], 5);
        resize(newData_rep2.setPos[s], 5);  
    }

    newData_rep1.setObs[0][0] = obsShort;
    newData_rep1.setPos[0][0] = 100;
    newData_rep1.setObs[0][1] = obsShort;
    newData_rep1.setPos[0][1] = 200;
    newData_rep1.setObs[0][2] = obsShort;
    newData_rep1.setPos[0][2] = 300;
    newData_rep1.setObs[0][3] = obsLong;
    newData_rep1.setPos[0][3] = 500;
    newData_rep1.setObs[0][4] = obsShort;
    newData_rep1.setPos[0][4] = 600;

    newData_rep2.setObs[0][0] = obsShort;
    newData_rep2.setPos[0][0] = 105;
    newData_rep2.setObs[0][1] = obsLong;
    newData_rep2.setPos[0][1] = 295;
    newData_rep2.setObs[0][2] = obsShort;
    newData_rep2.setPos[0][2] = 400;    
    newData_rep2.setObs[0][3] = obsShort;
    newData_rep2.setPos[0][3] = 505;
    newData_rep2.setObs[0][4] = obsShort;
    newData_rep2.setPos[0][4] = 595;    

    data_replicates[0] = newData_rep1;
    data_replicates[1] = newData_rep2;
}

void check(String<Data> &data_replicates)
{ 
    std::cout << "Rep1: \t";
    for (unsigned i = 0; i < length(data_replicates[0].setObs[0]); ++i)
    {
        std::cout << data_replicates[0].setPos[0][i] << " - " << (data_replicates[0].setPos[0][i] + data_replicates[0].setObs[0][i].length()) << "\t"; 
    }
    std::cout << std::endl;

    std::cout << "Rep2: \t";
    for (unsigned i = 0; i < length(data_replicates[1].setObs[0]); ++i)
    {
        std::cout << data_replicates[1].setPos[0][i] << " - " << (data_replicates[1].setPos[0][i] + data_replicates[1].setObs[0][i].length()) << "\t"; 
    }
    std::cout << std::endl;
}

void test_replicateIntersection()
{
    String<Data> data_replicates;
    setUp(data_replicates);

    check(data_replicates);
    
    intersect_replicateIntervals(data_replicates);

    check(data_replicates);
}
// END TEST
///////////////////////////////////////////// 

#endif
