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

template<typename TDOUBLE>
void prior_mle(GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2,  
               Data &data, 
               AppOptions &options)
{
    String<String<String<TDOUBLE> > > statePosteriors1;
    String<String<String<TDOUBLE> > > statePosteriors2;
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


template<typename TDOUBLE>
void prior_mle(GAMMA2_REG<TDOUBLE> &gamma1_reg, GAMMA2_REG<TDOUBLE> &gamma2_reg,  
               Data &data, 
               AppOptions &options)
{
    // compound gamma distribution: estimate parameters to use for initialization
    GAMMA2_REG<TDOUBLE> gammaC_reg(options.useKdeThreshold);
    gammaC_reg.k = 4.0;
    gammaC_reg.b0 = -1.5;   // TODO auto estimate !
    gammaC_reg.b1 = 0.9;

    // assign all sites with 1.0 to compound distribution
    String<String<String<TDOUBLE> > > statePosteriors;
    resize(statePosteriors, 2, Exact());
    // split into non-enriched and enriched
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(statePosteriors[s], length(data.setObs[s]), Exact());
        for(unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            resize(statePosteriors[s][i], data.setObs[s][i].length(), Exact());
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)
            {
                statePosteriors[s][i][t] = 1.0;
            }
        }
    }
    std::cout << "Estimate parameters of compound gamma distribution:" << std::endl;    
    // estimate parameters using GSL simplex2
    if (!gammaC_reg.updateRegCoeffsAndK(statePosteriors, data.setObs, options.g1_kMin, options.g1_kMax, options))
    {
        std::cerr << "ERROR: in updating parameters for compound gamma distribution using GSL simplex2." << std::endl;
    }
    myPrint(gammaC_reg);

    ///////////////////
    // classify sites as "non-enriched" or "enriched" based on regression parameter b1

    String<String<String<TDOUBLE> > > statePosteriors1;
    String<String<String<TDOUBLE> > > statePosteriors2;
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
                double x = std::max(data.setObs[s][i].rpkms[t], options.minRPKMtoFit);
                double pred = exp(gammaC_reg.b0 + gammaC_reg.b1 * x);

                if (data.setObs[s][i].kdes[t] < pred)                       // most probably non-enriched    .. (Note: pred + 2.0 for RNA-seq data)
                {                                                           // allow for variation! 
                    statePosteriors1[s][i][t] = 0.9;
                    statePosteriors2[s][i][t] = 0.1;
                }
                else                                                        // most probably enriched
                {
                    statePosteriors1[s][i][t] = 0.1;
                    statePosteriors2[s][i][t] = 0.9;
                }
            }
        }
    }
   
    gamma1_reg.k = 5.0;    //std::min(gammaC_reg.k + 5.0, options.g1_kMax); TODO why ??? 
    gamma2_reg.k = gammaC_reg.k;  
    gamma1_reg.b0 = gammaC_reg.b0 - 0.25;          // TODO ?
    gamma2_reg.b0 = gammaC_reg.b0 + 0.25; 
    gamma1_reg.b1 = gammaC_reg.b1; // + 0.1; 
    gamma2_reg.b1 = gammaC_reg.b1; // - 0.1; 

    std::cout << "Set parameters of gamma1 distribution:" << std::endl;    
    myPrint(gamma1_reg);
    std::cout << "Set parameters of gamma2 distribution:" << std::endl;    
    myPrint(gamma2_reg);
    
    // estimate parameters using GSL simplex2
    if (!gamma1_reg.updateRegCoeffsAndK(statePosteriors1, data.setObs, options.g1_kMin, options.g1_kMax, options))
    {
        std::cerr << "ERROR: in updating parameters for gamma1 distribution using GSL simplex2." << std::endl;
    }
    std::cout << "Updated parameters of gamma1 distribution using  GSL simplex2:" << std::endl;
    myPrint(gamma1_reg);
    
    if (!gamma2_reg.updateRegCoeffsAndK(statePosteriors2, data.setObs, options.g2_kMin, options.g2_kMax, options))
    {
        std::cerr << "ERROR: in updating parameters for gamma2 distribution using GSL simplex2." << std::endl;
    }
    std::cout << "Updated parameters of gamma2 distribution using  GSL simplex2:" << std::endl;
    myPrint(gamma2_reg);


    // TODO if user parameter given, set parameter to ...
    /*gamma1_reg.b0 = options.init_g1b0;   
    gamma2_reg.b0 = options.init_g2b0;
    gamma1_reg.b1 = options.init_g1b1; 
    gamma2_reg.b1 = options.init_g2b1;

    gamma1_reg.k = gamma1.k;
    gamma1_reg.tp = gamma1.tp;
    gamma2_reg.k = gamma2.k;
    gamma2_reg.tp = gamma2.tp;

    myPrint(gamma1_reg);
    myPrint(gamma2_reg);*/
}


template<typename TDOUBLE, typename TBIN>
void estimateTransitions(String<String<double> > &initTrans, 
                         GAMMA2<TDOUBLE> &gamma1, GAMMA2<TDOUBLE> &gamma2, TBIN &bin1, TBIN &bin2, 
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
            bool prev_valid = true;

            long double g1_d = 1.0;
            long double g2_d = 0.0;
            if (kde >= gamma1.tp) 
            {
                g1_d = gamma1.getDensity(kde);
                g2_d = gamma2.getDensity(kde); 
            }
            long double bin1_d = 1.0;
            long double bin2_d = 0.0;
            if (k > 0)
            {
                bin1_d = bin1.getDensity(k, n);
                bin2_d = bin2.getDensity(k, n);
            }

            double max = g1_d * bin1_d;   // most likely non-enriched and no crosslink "0"
            if (g1_d * bin2_d > max)      // most likely non-enriched and crosslink "1" 
            {
                max = g1_d * bin2_d;
                prev_state = 1;
            }
            if (g2_d * bin1_d > max)        // most likely enriched and no crosslink "2"
            {
                max = g2_d * bin1_d;
                prev_state = 2;
            }
            if (g2_d * bin2_d > max)        // most likely enriched and crosslink "3"
            {
                max = g2_d * bin2_d;
                prev_state = 3;
            }

            // count transitions
            for (unsigned t = 1; t < data.setObs[s][i].length(); ++t)
            {
                kde = data.setObs[s][i].kdes[t];
                k = data.setObs[s][i].truncCounts[t];
                n = data.setObs[s][i].nEstimates[t];
                unsigned curr_state = 0;

                g1_d = 1.0;
                g2_d = 0.0;
                if (kde >= gamma1.tp) 
                {
                    g1_d = gamma1.getDensity(kde);
                    g2_d = gamma2.getDensity(kde); 
                }
                bin1_d = 1.0;
                bin2_d = 0.0;
                if (k > 0)
                {
                    bin1_d = bin1.getDensity(k, n);
                    bin2_d = bin2.getDensity(k, n);
                }

                max = g1_d * bin1_d;          // most likely non-enriched and no crosslink "0"
                if (g1_d * bin2_d > max)      // most likely non-enriched and crosslink "1" 
                {
                    max = g1_d * bin2_d;
                    curr_state = 1;
                }
                if (g2_d * bin1_d > max)        // most likely enriched and no crosslink "2"
                {
                    max = g2_d * bin1_d;
                    curr_state = 2;
                }
                if (g2_d * bin2_d > max)        // most likely enriched and crosslink "3"
                {
                    max = g2_d * bin2_d;
                    curr_state = 3;
                }

                if (max > 0.0 && prev_valid) 
                    ++transFreqs[prev_state][curr_state];
                //else
                //    std::cout << "Note: estimating transFreqs, all emission probs 0!" << gamma1.getDensity(kde) << " - " << gamma2.getDensity(kde) << " - " <<  bin1_d<< " - " << bin2_d<< "  n: " << n << " k: " << k << std::endl;
                
                if (max > 0.0) 
                    prev_valid = true;
                else
                    prev_valid = false;

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


template<typename TDOUBLE, typename TBIN>
void estimateTransitions(String<String<double> > &initTrans, 
                         GAMMA2_REG<TDOUBLE> &gamma1, GAMMA2_REG<TDOUBLE> &gamma2, TBIN &bin1, TBIN &bin2, 
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
            double d1_pred = exp(gamma1.b0 + gamma1.b1 * data.setObs[s][i].rpkms[0]);
            double d2_pred = exp(gamma2.b0 + gamma2.b1 * data.setObs[s][i].rpkms[0]);
            unsigned prev_state = 0;
            bool prev_valid = true;

            long double g1_d = 1.0;
            long double g2_d = 0.0;
            if (kde >= gamma1.tp)
            {
                g1_d = gamma1.getDensity(kde, d1_pred, options);
                g2_d = gamma2.getDensity(kde, d2_pred, options);
            }
            long double bin1_d = 1.0;
            long double bin2_d = 0.0;
            if (k > 0)
            {
                bin1_d = bin1.getDensity(k, n);
                bin2_d = bin2.getDensity(k, n);
            }

            long double max = g1_d * bin1_d;   // most likely non-enriched and no crosslink "0"
            if (g1_d * bin2_d> max)      // most likely non-enriched and crosslink "1" 
            {
                max = g1_d * bin2_d;
                prev_state = 1;
            }
            if (g2_d * bin1_d> max)        // most likely enriched and no crosslink "2"
            {
                max = g2_d * bin1_d;
                prev_state = 2;
            }
            if (g2_d * bin2_d> max)        // most likely enriched and crosslink "3"
            {
                max = g2_d * bin2_d;
                prev_state = 3;
            }
            // count transitions
            for (unsigned t = 1; t < data.setObs[s][i].length(); ++t)
            {
                kde = data.setObs[s][i].kdes[t];
                k = data.setObs[s][i].truncCounts[t];
                n = data.setObs[s][i].nEstimates[t];
                d1_pred = exp(gamma1.b0 + gamma1.b1 * data.setObs[s][i].rpkms[t]);
                d2_pred = exp(gamma2.b0 + gamma2.b1 * data.setObs[s][i].rpkms[t]);
                unsigned curr_state = 0;

                g1_d = 1.0;
                g2_d = 0.0;
                if (kde >= gamma1.tp)
                {
                    g1_d = gamma1.getDensity(kde, d1_pred, options);
                    g2_d = gamma2.getDensity(kde, d2_pred, options);
                }
                bin1_d = 1.0;
                bin2_d = 0.0;
                if (k > 0)
                {
                    bin1_d = bin1.getDensity(k, n);
                    bin2_d = bin2.getDensity(k, n);
                }

                max = g1_d * bin1_d;          // most likely non-enriched and no crosslink "0"
                if (g1_d * bin2_d> max)      // most likely non-enriched and crosslink "1" 
                {
                    max = g1_d * bin2_d;
                    curr_state = 1;
                }
                if (g2_d * bin1_d> max)        // most likely enriched and no crosslink "2"
                {
                    max = g2_d * bin1_d;
                    curr_state = 2;
                }
                if (g2_d * bin2_d> max)        // most likely enriched and crosslink "3"
                {
                    max = g2_d * bin2_d;
                    curr_state = 3;
                }
                if (max > 0.0 && prev_valid) 
                    ++transFreqs[prev_state][curr_state];
                //else
                //    std::cout << "Note: estimating transFreqs, all emission probs 0!" << gamma1.getDensity(kde) << " - " << gamma2.getDensity(kde) << " - " <<  bin1_d<< " - " << bin2_d<< std::endl;
                
                if (max > 0.0) 
                    prev_valid = true;
                else
                    prev_valid = false;

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
