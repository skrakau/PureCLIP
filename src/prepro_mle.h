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

#ifndef APPS_HMMS_PREPRO_MLE_H_
#define APPS_HMMS_PREPRO_MLE_H_
  
#include <iostream>
#include <fstream>

using namespace seqan;


void prior_mle(GAMMA2 &gamma1, GAMMA2 &gamma2,  
               Data &data, 
               AppOptions &options)
{
    String<String<String<long double> > > statePosteriors1;
    String<String<String<long double> > > statePosteriors2;
    resize(statePosteriors1, 2, Exact());
    resize(statePosteriors2, 2, Exact());
    // split into non-enriched and enriched
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors1[s], length(data.setObs[s]), Exact());
        resize(statePosteriors2[s], length(data.setObs[s]), Exact());
        for(unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            resize(statePosteriors1[s][i], data.setObs[s][i].length(), Exact());
            resize(statePosteriors2[s][i], data.setObs[s][i].length(), Exact());
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)
            {
                if (data.setObs[s][i].kdes[t] < options.prior_kdeThreshold)    // most probably non-enriched
                {                                                           // use 1.0; 0.0; to avoid bias by large number of zeros 
                    statePosteriors1[s][i][t] = 0.999;
                    statePosteriors2[s][i][t] = 0.001;
                }
                else                                                        // most probably enriched
                {
                    statePosteriors1[s][i][t] = 0.001;
                    statePosteriors2[s][i][t] = 0.999;
                }
            }
        }
    }

    gamma1.k = 0.50; 
    gamma1.updateTheta(statePosteriors1, data.setObs, options);
    gamma1.updateK(statePosteriors1, data.setObs, options.g1_kMin, options.g1_kMax, options);
    gamma2.k = 1.5; 
    gamma2.updateTheta(statePosteriors2, data.setObs, options);
    gamma2.updateK(statePosteriors2, data.setObs, options.g2_kMin, options.g2_kMax, options);

    if (options.verbosity >= 1) 
    {
        std::cout << "Initialization:" << std::endl;
        std::cout << "Gamma1 theta: " << gamma1.theta << std::endl;
        std::cout << "Gamma1 k: " << gamma1.k << std::endl;
        std::cout << "Gamma2 theta: " << gamma2.theta << std::endl;
        std::cout << "Gamma2 k: " << gamma2.k << std::endl;
    }

    myPrint(gamma1);
    myPrint(gamma2);
}

void prior_mle(GAMMA2_REG &gamma1_reg, GAMMA2_REG &gamma2_reg,  
               Data &data, 
               AppOptions &options)
{
    String<String<String<long double> > > statePosteriors1;
    String<String<String<long double> > > statePosteriors2;
    resize(statePosteriors1, 2, Exact());
    resize(statePosteriors2, 2, Exact());
    // split into non-enriched and enriched
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors1[s], length(data.setObs[s]), Exact());
        resize(statePosteriors2[s], length(data.setObs[s]), Exact());
        for(unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            resize(statePosteriors1[s][i], data.setObs[s][i].length(), Exact());
            resize(statePosteriors2[s][i], data.setObs[s][i].length(), Exact());
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)
            {
                if (data.setObs[s][i].kdes[t] < options.prior_kdeThreshold)    // most probably non-enriched
                {                                                           // use 1.0; 0.0; to avoid bias by large number of zeros 
                    statePosteriors1[s][i][t] = 0.999;
                    statePosteriors2[s][i][t] = 0.001;
                }
                else                                                        // most probably enriched
                {
                    statePosteriors1[s][i][t] = 0.001;
                    statePosteriors2[s][i][t] = 0.999;
                }
            }
        }
    }

    GAMMA2     gamma1(options.useKdeThreshold);  
    GAMMA2     gamma2(options.useKdeThreshold);      
    gamma1.k = 1.0; 
    gamma1.updateTheta(statePosteriors1, data.setObs, options);
    gamma1.updateK(statePosteriors1, data.setObs, options.g1_kMin, options.g1_kMax,  options);
    gamma2.k = 1.5;  
    gamma2.updateTheta(statePosteriors2, data.setObs, options);
    gamma2.updateK(statePosteriors2, data.setObs, options.g2_kMin, options.g2_kMax, options);

    if (options.verbosity >= 1) 
    {
        std::cout << "Initialization:" << std::endl;
        std::cout << "Gamma1 theta: " << gamma1.theta << std::endl;
        std::cout << "Gamma1 k: " << gamma1.k << std::endl;
        std::cout << "Gamma2 theta: " << gamma2.theta << std::endl;
        std::cout << "Gamma2 k: " << gamma2.k << std::endl;
    }

    myPrint(gamma1);
    myPrint(gamma2);

    gamma1_reg.mean = gamma1.theta * gamma1.k;
    gamma2_reg.mean = gamma2.theta * gamma2.k;

    gamma1_reg.b0 = log(gamma1_reg.mean);
    gamma2_reg.b0 = log(gamma2_reg.mean);


    gamma1_reg.k = gamma1.k;
    gamma1_reg.tp = gamma1.tp;
    gamma2_reg.k = gamma2.k;
    gamma2_reg.tp = gamma2.tp;

    myPrint(gamma1_reg);
    myPrint(gamma2_reg);
}




template<typename TB1, typename TB2>
void estimateTransitions(String<String<double> > &initTrans, 
                         GAMMA2 &gamma1, GAMMA2 &gamma2, TB1 &bin1, TB2 &bin2, 
                         Data &data,
                         AppOptions &options)
{
    // split into non-enriched and enriched
    String<String<unsigned> >  transFreqs;
    resize(transFreqs, 4, Exact());
    resize(initTrans, 4, Exact());
    for (unsigned k = 0; k < 4; ++k)
    {
        resize(transFreqs[k], 4, 0, Exact());
        resize(initTrans[k], 4, Exact());
    }

    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)  
        {
            double kde = data.setObs[s][i].kdes[0];
            unsigned k = data.setObs[s][i].truncCounts[0];
            unsigned n = data.setObs[s][i].nEstimates[0];
            unsigned prev_state = 0;
            double max = gamma1.getDensity(kde) * bin1.getDensity(k, n);   // most likely non-enriched and no crosslink "0"
         
            if (gamma1.getDensity(kde) * bin2.getDensity(k, n) > max)      // most likely non-enriched and crosslink "1" 
            {
                max = gamma1.getDensity(kde) * bin2.getDensity(k, n);
                prev_state = 1;
            }
            if (gamma2.getDensity(kde) * bin1.getDensity(k, n) > max)        // most likely enriched and no crosslink "2"
            {
                max = gamma2.getDensity(kde) * bin1.getDensity(k, n);
                prev_state = 2;
            }
            if (gamma2.getDensity(kde) * bin2.getDensity(k, n) > max)        // most likely enriched and crosslink "3"
            {
                max = gamma2.getDensity(kde) * bin2.getDensity(k, n);
                prev_state = 3;
            }

            // count transitions
            for (unsigned t = 1; t < data.setObs[s][i].length(); ++t)
            {
                kde = data.setObs[s][i].kdes[t];
                k = data.setObs[s][i].truncCounts[t];
                n = data.setObs[s][i].nEstimates[t];
                unsigned curr_state = 0;
                max = gamma1.getDensity(kde) * bin1.getDensity(k, n);          // most likely non-enriched and no crosslink "0"

                if (gamma1.getDensity(kde) * bin2.getDensity(k, n) > max)      // most likely non-enriched and crosslink "1" 
                {
                    max = gamma1.getDensity(kde) * bin2.getDensity(k, n);
                    curr_state = 1;
                }
                if (gamma2.getDensity(kde) * bin1.getDensity(k, n) > max)        // most likely enriched and no crosslink "2"
                {
                    max = gamma2.getDensity(kde) * bin1.getDensity(k, n);
                    curr_state = 2;
                }
                if (gamma2.getDensity(kde) * bin2.getDensity(k, n) > max)        // most likely enriched and crosslink "3"
                {
                    max = gamma2.getDensity(kde) * bin2.getDensity(k, n);
                    curr_state = 3;
                }
                ++transFreqs[prev_state][curr_state];
                prev_state = curr_state;
            }
        }
    }
    // convert counts to transition probabilities
    for (unsigned k_1 = 0; k_1 < 4; ++k_1)
    {
        unsigned sum = 0;   // pseudocount
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
            sum += ++transFreqs[k_1][k_2];
   
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
        {
            initTrans[k_1][k_2] = (double)transFreqs[k_1][k_2] / (double)sum;     
        }
    }
    if (options.verbosity >= 1) 
    {
        std::cout << " initial transition frequencies: " << std::endl;
        for (unsigned k_1 = 0; k_1 < 4; ++k_1)
        {
            std::cout << "  " << k_1 << ": ";
            for (unsigned k_2 = 0; k_2 < 4; ++k_2)
                std::cout << transFreqs[k_1][k_2] << "  ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}


template<typename TB1, typename TB2>
void estimateTransitions(String<String<double> > &initTrans, 
                         GAMMA2_REG &gamma1, GAMMA2_REG &gamma2, TB1 &bin1, TB2 &bin2, 
                         Data &data,
                         AppOptions &options)
{
    // split into non-enriched and enriched
    String<String<unsigned> >  transFreqs;
    resize(transFreqs, 4, Exact());
    resize(initTrans, 4, Exact());
    for (unsigned k = 0; k < 4; ++k)
    {
        resize(transFreqs[k], 4, 0, Exact());
        resize(initTrans[k], 4, Exact());
    }

    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)  
        {
            double d1_pred = exp(gamma1.b0 + gamma1.b1 * data.setObs[s][i].rpkms[0]);
            double d2_pred = exp(gamma2.b0 + gamma2.b1 * data.setObs[s][i].rpkms[0]);
            unsigned prev_state = 0;
            long double gamma1_eProb = gamma1.getDensity(data.setObs[s][i].kdes[0], d1_pred);
            long double gamma2_eProb = gamma2.getDensity(data.setObs[s][i].kdes[0], d2_pred);
            if (data.setObs[s][i].kdes[0] < gamma1.tp)
            {
                gamma1_eProb = 1.0;
                gamma2_eProb = 0.0; 
            }
            unsigned k = data.setObs[s][i].truncCounts[0];
            unsigned n = data.setObs[s][i].nEstimates[0];

            long double max = gamma1_eProb * bin1.getDensity(k, n);   // most likely non-enriched and no crosslink "0"
         
            if (gamma1_eProb * bin2.getDensity(k, n) > max)      // most likely non-enriched and crosslink "1" 
            {
                max = gamma1_eProb * bin2.getDensity(k, n);
                prev_state = 1;
            }
            if (gamma2_eProb * bin1.getDensity(k, n) > max)        // most likely enriched and no crosslink "2"
            {
                max = gamma2_eProb * bin1.getDensity(k, n);
                prev_state = 2;
            }
            if (gamma2_eProb * bin2.getDensity(k, n) > max)        // most likely enriched and crosslink "3"
            {
                max = gamma2_eProb * bin2.getDensity(k, n);
                prev_state = 3;
            }
            // count transitions
            for (unsigned t = 1; t < data.setObs[s][i].length(); ++t)
            {
                d1_pred = exp(gamma1.b0 + gamma1.b1 * data.setObs[s][i].rpkms[t]);
                d2_pred = exp(gamma2.b0 + gamma2.b1 * data.setObs[s][i].rpkms[t]);
                k = data.setObs[s][i].truncCounts[t];
                n = data.setObs[s][i].nEstimates[t];

                unsigned curr_state = 0;
                gamma1_eProb = gamma1.getDensity(data.setObs[s][i].kdes[t], d1_pred);
                gamma2_eProb = gamma2.getDensity(data.setObs[s][i].kdes[t], d2_pred);
                if (data.setObs[s][i].kdes[t] < gamma1.tp)
                {
                    gamma1_eProb = 1.0;
                    gamma2_eProb = 0.0; 
                }

                max = gamma1_eProb * bin1.getDensity(k, n);          // most likely non-enriched and no crosslink "0"

                if (gamma1_eProb * bin2.getDensity(k, n) > max)      // most likely non-enriched and crosslink "1" 
                {
                    max = gamma1_eProb * bin2.getDensity(k, n);
                    curr_state = 1;
                }
                if (gamma2_eProb * bin1.getDensity(k, n) > max)        // most likely enriched and no crosslink "2"
                {
                    max = gamma2_eProb * bin1.getDensity(k, n);
                    curr_state = 2;
                }
                if (gamma2_eProb * bin2.getDensity(k, n) > max)        // most likely enriched and crosslink "3"
                {
                    max = gamma2_eProb * bin2.getDensity(k, n);
                    curr_state = 3;
                }
                ++transFreqs[prev_state][curr_state];
                prev_state = curr_state;
            }
        }
    }
    // convert counts to transition probabilities
    for (unsigned k_1 = 0; k_1 < 4; ++k_1)
    {
        unsigned sum = 0;   // pseudocount
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
            sum += ++transFreqs[k_1][k_2];
   
        for (unsigned k_2 = 0; k_2 < 4; ++k_2)
        {
            initTrans[k_1][k_2] = (double)transFreqs[k_1][k_2] / (double)sum;     
        }
    }
 
    if (options.verbosity >= 1) 
    {
        std::cout << " initial transition frequencies: " << std::endl;
        for (unsigned k_1 = 0; k_1 < 4; ++k_1)
        {
            std::cout << "  " << k_1 << ": ";
            for (unsigned k_2 = 0; k_2 < 4; ++k_2)
                std::cout << transFreqs[k_1][k_2] << "  ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}



    
#endif
