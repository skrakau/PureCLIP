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



#ifndef APPS_HMMS_CALL_SITES_H_
#define APPS_HMMS_CALL_SITES_H_

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
#include "call_sites_replicates.h"

#include "density_functions_reg.h"
#include "density_functions_crosslink.h"
#include "density_functions_crosslink_reg.h"

using namespace seqan;

// intervals: extract intervals from input string and get store contigIds 
template <typename TStore>
bool loadIntervals(AppOptions &options, TStore &store)
{
    typedef StringSet<CharString>   TNameStore;

    NameStoreCache<TNameStore>  nameStoreCache(store.contigNameStore);
    if (options.verbosity >= 1 && !empty(options.intervals_str)) 
        std::cout << "Load user specified genomic intervals for learning ..." << std::endl;
    for (unsigned i = 0; i < length(options.intervals_str); ++i)
    {
        CharString buffer;
        while (i < length(options.intervals_str) && options.intervals_str[i] != ';')
        {
            appendValue(buffer, options.intervals_str[i]);
            ++i;
        }
        unsigned j = 0;
        CharString contigName;
        while (j < length(buffer) && buffer[j] != ':')
        {
            appendValue(contigName, buffer[j]);
            ++j;
            //++i;
        }
        ++j;
        //++i;
        // retrieve corresponding contigId
        unsigned contigId;
        if (!getIdByName(contigId, nameStoreCache, contigName))
        {
            std::cerr << "ERROR: Can't find contigName in store.contigNameStore!" << contigName << std::endl;
            continue;
        }

        appendValue(options.intervals_contigIds, contigId);
        // extract positions
        unsigned i1 = 0;
        unsigned i2 = length(store.contigStore[contigId].seq);
        if (j < length(buffer))
        {
            CharString i1_str;
            while (buffer[j] != '-')
            {
                appendValue(i1_str, buffer[j]);
                ++j;
            }
            SEQAN_ASSERT_EQ(buffer[j], '-');
            CharString i2_str;
            while (j < length(buffer))
            {
                appendValue(i2_str, buffer[j]);
                ++j;
            } 
            i1 = atoi(toCString(i1_str));
            i2 = atoi(toCString(i2_str));
        }
        if (options.verbosity >= 2) std::cout << "contigName: " << contigName << " i1: " << i1 << " i2: " << i2 << std::endl;

        String<unsigned> interval;
        appendValue(interval, i1);
        appendValue(interval, i2);
        appendValue(options.intervals_positions, interval);
    }

    if (empty(options.intervals_str))   // use all reference contigs
    {
        for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId)
        {
            appendValue(options.intervals_contigIds, contigId);
            CharString contigName = store.contigNameStore[contigId];

            unsigned i1 = 0;
            unsigned i2 = length(store.contigStore[contigId].seq);
            if (options.verbosity >= 2)  std::cout << "contigName: " << contigName << " i1: " << i1 << " i2: " << i2 << std::endl;

            String<unsigned> interval;
            appendValue(interval, i1);
            appendValue(interval, i2);
            appendValue(options.intervals_positions, interval);
        }
    }

    return 0;
}


template <typename TStore>
bool loadApplyChrs(AppOptions &options, TStore &store)
{
    typedef StringSet<CharString>   TNameStore;

    NameStoreCache<TNameStore>  nameStoreCache(store.contigNameStore);
    if (options.verbosity >= 1 && !empty(options.applyChr_str)) 
        std::cout << "Load user specified chromosomes for applying HMM model ..." << std::endl;
    for (unsigned i = 0; i < length(options.applyChr_str); ++i)
    {
        CharString contigName;
        while (i < length(options.applyChr_str) && options.applyChr_str[i] != ';')
        {
            appendValue(contigName, options.applyChr_str[i]);
            ++i;
        }
        // retrieve corresponding contigId
        unsigned contigId;
        if (!getIdByName(contigId, nameStoreCache, contigName))
        {
            std::cerr << "ERROR: Can't find contigName in store.contigNameStore!" << contigName << std::endl;
            continue;
        }
        if (options.verbosity >= 2) std::cout << "contigName: " << contigName << std::endl;
        appendValue(options.applyChr_contigIds, contigId);
    }

    if (empty(options.applyChr_str))   // use all reference contigs
    {
        for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId)
        {
            if (options.verbosity >= 2) std::cout << "contigName: " << store.contigNameStore[contigId] << std::endl;
            appendValue(options.applyChr_contigIds, contigId);
        }
    }

    return 0;
}



template <typename TContigObservations, typename TBai, typename TStore>
int loadObservations(TContigObservations &contigObservationsF, TContigObservations &contigObservationsR, unsigned contigId, CharString bamFileName, TBai &baiIndex, TStore &store, AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    if (options.verbosity >= 2) std::cout << "Parse alignments ... " << std::endl;

    // Open BamFileIn for reading.
    if (options.verbosity >= 2) std::cout << "Open Bam and Bai file ... "  << "\n";
    BamFileIn inFile;
    if (!open(inFile, toCString(bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << bamFileName << " for reading.\n";
        return 1;
    }
    BamHeader header;
    readHeader(header, inFile);

    // Translate from contig name to rID.
    int rID = 0;
    if (!getIdByName(rID, contigNamesCache(context(inFile)), store.contigNameStore[contigId]))
    {
        if (options.verbosity >= 2) std::cout << "NOTE: Contig " << store.contigNameStore[contigId] << " not existing in BAM file.\n";
        return 2; 
    }

    resize(contigObservationsF.truncCounts, length(store.contigStore[contigId].seq), 0, Exact());
    resize(contigObservationsR.truncCounts, length(store.contigStore[contigId].seq), 0, Exact());

    if (!parse_bamRegion(contigObservationsF, contigObservationsR, inFile, baiIndex, rID, options))
        return 2;

    // ATTENTIONE: reverse in-place here to avoid problems for observations datastructures (and use Modifier iterator later within writeStates)  !!!!!!!!
    reverse(contigObservationsR);      

    if (options.verbosity >= 2) std::cout << "... observations loaded" << std::endl;
#ifdef HMM_PROFILE
    Times::instance().time_loadObservations += (sysTime() - timeStamp);
#endif
    return 0;
}


template <typename TBai, typename TStore, typename TOptions>
bool loadBAMCovariates(Data &data, TBai &inputBaiIndex, TStore &store, bool parallelize, TOptions &options)
{
    if (options.verbosity >= 1) std::cout << "  Parse input BAM, get truncCounts, compute KDEs ... " << std::endl;

    bool stop = false;
    // TODO check if working
#if HMM_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(parallelize ? 2 : 1))
#endif 
    for (unsigned s = 0; s < 2; ++s)
    {
        // Open BamFileIn for reading.
        SEQAN_OMP_PRAGMA(critical)
        std::cout << "Open Bam and Bai file ... (strand: " << s << ")\n";
        BamFileIn inFile;
        if (!open(inFile, toCString(options.inputBamFileName)))
        {
            SEQAN_OMP_PRAGMA(critical)
            std::cerr << "ERROR: Could not open " << options.inputBamFileName << " for reading.\n";
            stop = true;
            continue;
        }
        BamHeader header;
        readHeader(header, inFile);

        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            // Translate from contig name to rID.
            int rID = 0;
            if (!getIdByName(rID, contigNamesCache(context(inFile)), store.contigNameStore[data.setObs[s][i].contigId]))
            {
                SEQAN_OMP_PRAGMA(critical)
                std::cerr << "ERROR: Contig " << store.contigNameStore[data.setObs[s][i].contigId] << " not known.\n";
                stop = true;
                continue;
            }
            String<__uint16> truncCounts;
            resize(truncCounts, data.setObs[s][i].length(), 0, Exact());
            if (s == 0)
            {
                unsigned beginPos = data.setPos[s][i];
                unsigned endPos = data.setPos[s][i] + data.setObs[s][i].length();
                parse_bamRegion(truncCounts, inFile, inputBaiIndex, rID, beginPos, endPos, true, options);
            }
            else
            {
                unsigned beginPos = length(store.contigStore[data.setObs[s][i].contigId].seq) - (data.setObs[s][i].length() + data.setPos[s][i]);
                unsigned endPos = beginPos + data.setObs[s][i].length();

                parse_bamRegion(truncCounts, inFile, inputBaiIndex, rID, beginPos, endPos, false, options);
                // reverse
                reverse(truncCounts);
            }
            // compute KDEs
            data.setObs[s][i].computeKDEs(truncCounts, options);
        }
    }
    if (stop) return false;

    if (options.verbosity >= 2) std::cout << "... BAM covariates loaded" << std::endl;
    return true;
}


// precomputed 
template <typename TStore>
bool loadCovariates(String<double> &contigCovsF, String<double> &contigCovsR, unsigned contigId, TStore &store, AppOptions &options)
{   
    if (!empty(options.rpkmFileName)) 
    {
        double minRPKM = 0.0;
        if (options.useLogRPKM) 
            minRPKM = options.minRPKMtoFit - 1.0; 
        
        resize(contigCovsF, length(store.contigStore[contigId].seq), minRPKM, Exact());
        resize(contigCovsR, length(store.contigStore[contigId].seq), minRPKM, Exact());

        if (options.verbosity >= 1) std::cout << "Parse covariates ... input signal" << std::endl;
        BedFileIn bedIn(toCString(options.rpkmFileName));
        BedRecord<Bed6> bedRecord;      
        while (!atEnd(bedIn))
        {
            try
            {
                readRecord(bedRecord, bedIn);
            }
            catch (ParseError const & e)
            {
                std::cerr << "ERROR: input BED record is badly formatted. " << e.what() << "\n";
            }
            catch (IOError const & e)
            {
                std::cerr << "ERROR: could not copy input BED record. " << e.what() << "\n";
            } 

            if (bedRecord.ref == store.contigNameStore[contigId])   
            {
               if (bedRecord.strand == '+')
               {
                    std::stringstream ss(toCString(bedRecord.score)); 
                    double score;
                    ss >> score;
                    if (options.useLogRPKM  && score > 0.0)
                        score = log(score);
                    else if (options.useLogRPKM)
                        score = minRPKM;

                    for (int i = bedRecord.beginPos; i < bedRecord.endPos; ++i)
                    {
                        if (contigCovsF[i] == minRPKM)
                            contigCovsF[i] = score;
                        else
                            contigCovsF[i] = std::max(score, contigCovsF[i]);
                    }
               }
               else
               {
                    std::stringstream ss(toCString(bedRecord.score)); 
                    double score;
                    ss >> score;
                    if (options.useLogRPKM && score > 0.0)
                        score = log(score);
                    else if (options.useLogRPKM)
                        score = minRPKM;

                    for (int i = bedRecord.beginPos; i < bedRecord.endPos; ++i)
                    {
                         if (contigCovsR[i] == minRPKM)
                            contigCovsR[i] = score;
                        else
                            contigCovsR[i] = std::max(score, contigCovsR[i]);
                   }
               }
            }
        } 

        // ATTENTIONE: reverse in-place here to avoid problems for observations datastructures (and use Modifier iterator later within writeStates)  !!!!!!!!
        reverse(contigCovsR);      

        if (options.verbosity >= 2) std::cout << "... covariates loaded" << std::endl;
    }
    return 0;
}


template <typename TStore>
bool loadMotifCovariates(String<String<float> > &contigCovs, String<String<char> > &motifIds, unsigned contigId, TStore &store, AppOptions &options)
{   
    
    resize(contigCovs, 2, Exact());
    resize(motifIds, 2, Exact());
    for (unsigned s = 0; s < 2; ++s)
    {
        resize(contigCovs[s], length(store.contigStore[contigId].seq), 0.0, Exact());
        resize(motifIds[s], length(store.contigStore[contigId].seq), 0, Exact());
    }

    if (!empty(options.fimoFileName)) 
    {
        if (options.verbosity >= 2) std::cout << "Parse covariates ... fimo scores for input motifs" << std::endl;
        BedFileIn bedIn(toCString(options.fimoFileName));
        BedRecord<Bed6> bedRecord;      
        while (!atEnd(bedIn))
        {
            try
            {
                readRecord(bedRecord, bedIn);
            }
            catch (ParseError const & e)
            {
                std::cerr << "ERROR: input BED record is badly formatted. " << e.what() << "\n";
            }
            catch (IOError const & e)
            {
                std::cerr << "ERROR: could not copy input BED record. " << e.what() << "\n";
            } 

            if (bedRecord.ref == store.contigNameStore[contigId])   
            {
                std::stringstream ss(toCString(bedRecord.score)); 
                double score;
                ss >> score;

                // assume only one motif with one score for each position!
                // assume id 1-based
                unsigned id = atoi(toCString(bedRecord.name)) - 1;     
                if (id < options.nInputMotifs)
                {
                    if (bedRecord.strand == '+')
                    {
                        if (bedRecord.beginPos < (int)length(contigCovs[0]))
                        {
                            contigCovs[0][bedRecord.beginPos] = std::max(score, 0.0);         // ignore negative scores for the moment
                            motifIds[0][bedRecord.beginPos] = id;
                        }
                        else
                            std::cout << "Warning: beginPos of fimo score is not within contig! Ignored. (contigName: " << bedRecord.ref << ", beginPos: " << bedRecord.beginPos << ")" << std::endl;
                    }
                    else
                    {
                        if (bedRecord.beginPos < (int)length(contigCovs[1]))
                        {
                            contigCovs[1][bedRecord.beginPos] = std::max(score, 0.0);
                            motifIds[1][bedRecord.beginPos] = id;
                        }
                        else
                            std::cout << "Warning: beginPos of fimo score is not within contig! Ignored. (contigName: " << bedRecord.ref << ", beginPos: " << bedRecord.beginPos << ")" << std::endl;
                    }
                }
            }
        } 

        // ATTENTIONE: reverse in-place here to avoid problems for observations datastructures (and use Modifier iterator later within writeStates)  !!!!!!!!
        reverse(contigCovs[1]);    
        reverse(motifIds[1]);  

        if (options.verbosity >= 2) std::cout << "... fimo input motif score covriates loaded" << std::endl;
    }
    return 0;
}


template<typename TSeq, typename TOptions>
bool checkForPolyA(TSeq const& seq, TOptions &options)
{
     typedef typename Iterator<TSeq, Rooted>::Type TIter;

     TIter it1 = begin(seq);
     TIter it2 = begin(seq);    

     while (it1 != end(seq))
     {
         it2 = it1;
         while( it2 != end(seq) && (ordValue(*it2) == 0)) ++it2;

         if ((position(it2) - position(it1)) >= options.polyAThreshold) return true;           // use options value
         if (it2 != it1) 
             it1 = it2;
         else
             ++it1;
     }
     return false;
}


template<typename TSeq, typename TOptions>
bool checkForPolyT(TSeq const& seq, TOptions &options)
{
     typedef typename Iterator<TSeq, Rooted>::Type TIter;

     TIter it1 = begin(seq);
     TIter it2 = begin(seq);

     while (it1 != end(seq))
     {
         it2 = it1;
         while( it2 != end(seq) && (ordValue(*it2) == 3 )) ++it2;

         if ((position(it2) - position(it1)) >= options.polyAThreshold) return true;           // use options value
         if (it2 != it1) 
             it1 = it2;
         else
             ++it1;
     }
     return false;
}


template <typename TOptions>
void cleanCoveredIntervals(Data &data, unsigned contigLength, bool learning, TOptions &options)
{
    for (unsigned s = 0; s < 2; ++s)
    {
        String<Observations> tmp_setObs;
        String<unsigned> tmp_setPos;
        resize(tmp_setObs, length(data.setObs[s]));
        resize(tmp_setPos, length(data.setPos[s]));

        unsigned j = 0;
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            unsigned count = 0;
            bool discard = false;
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)
            {
                // do not consider intervals containing sites with extremely high read start counts at one position for learning 
                // - possible artifacts (PCR duplicates, mapping artifacts)
                // - causing problems when getting higher for likelihood computations
                // - indiv. counts limited to maxTruncCount2  
                if (learning && (data.setObs[s][i].truncCounts[t] >= options.maxTruncCount))
                {
                    discard = true;
                    if (options.verbosity >= 2)
                    {
                        if (s == 0)
                            std::cout << "Note: Ignore interval for learning due to position with >= " << options.maxTruncCount << " read starts at position " << (t + data.setPos[s][i]) << "." << std::endl;
                        else
                            std::cout << "Note: Ignore interval for learning due to position with >= " << options.maxTruncCount << " read starts at position " << (contigLength - (t + data.setPos[s][i]) - 1) << "." << std::endl;
                    }
                    break;
                }
                count += data.setObs[s][i].truncCounts[t];
            }

            if ((!options.discardSingletonIntervals || count > 1) && !discard)
            {
                tmp_setObs[j] = data.setObs[s][i];      
                tmp_setPos[j] = data.setPos[s][i];
                ++j;
            }
        }
        erase(tmp_setObs, j, length(tmp_setObs));
        erase(tmp_setPos, j, length(tmp_setPos));

        data.setObs[s] = tmp_setObs;
        data.setPos[s] = tmp_setPos;
    }
}


template <typename TContigObservations, typename TStore, typename TOptions>
void extractCoveredIntervals(Data &data, 
                             TContigObservations &contigObservationsF, TContigObservations &contigObservationsR,
                             String<double> &contigCovsF, String<double> &contigCovsR,
                             String<String<float> > &contigCovsFimo, 
                             String<String<char> > &motifIds, 
                             unsigned contigId, 
                             unsigned i1, unsigned i2,
                             bool excludePolyA,
                             bool excludePolyT,
                             bool learning,
                             TStore &store,
                             TOptions &options)
{
    // learn HMM for specified interval

    unsigned countPolyAs = 0;
    unsigned countPolyTs = 0;
    // FORWARD                                      // TODO merge code F and R!
    unsigned c1;                                   // covered interval begin
    unsigned c2;                                   // covered interval end
    unsigned i = i1;
    unsigned prev_c1 = i1;
    unsigned prev_c2 = i1;
    bool prev_dis = false;
    if (options.verbosity >= 2) std::cout << "F: Parse covered intervals and get observations  ..." << "i1: " << i1 << " i2: " << i2 <<  std::endl;
    while (i < i2 && (prev_c2 < i2 && !prev_dis))
    {

        while (i < i2 && contigObservationsF.truncCounts[i] == 0) ++i;    // find begin of covered interval     
        c1 = i;
        //std::cout << "TEST: i " << i << std::endl;
        if (((int)c1 - (int)options.intervalOffset) > (int)prev_c2)    // if gap bigger than intervalOffset, shift c1 to left
        {
            c1 -= options.intervalOffset;
        }
        else if (prev_c2 == i1)             // else if begin of user defined interval. set c1 to this
        {
            c1 = i1;
        }
        else if (prev_dis)                  // else if prev. interval was already discarded (and hence not appended)
        {
            c1 = prev_c1;
            ++i;
            if (i >= i2) break;
                
            while (i < i2 && contigObservationsF.truncCounts[i] > 0) ++i;      // find end of covered interval
            c2 = std::min(i + options.intervalOffset, i2);
            prev_c2 = c2;
            continue;
        }
        else       
        {                                   // else merge with previous interval 
            c1 = prev_c1;
            //std::cout << "TEST: eraseBack .. set c1 " << c1 << std::endl;
            eraseBack(data.setObs[0]); 
            eraseBack(data.setPos[0]);    
        }
        prev_c1 = c1;
        ++i;

        if (i >= i2) break;
            
        while (i < i2 && contigObservationsF.truncCounts[i] > 0) ++i;      // find end of covered interval
        c2 = std::min(i + options.intervalOffset, i2);

        if (excludePolyA) // check if covered interval contains internal polyA  
        {
            if (checkForPolyA(infix(store.contigStore[contigId].seq, c1, c2), options)) 
            {
                ++countPolyAs;
                prev_dis = true;
                prev_c2 = c2;
                continue;
            }
        }
        if (excludePolyT) // check for polyT (polyU) 
        {
            if (checkForPolyT(infix(store.contigStore[contigId].seq, c1, c2), options)) 
            {
                ++countPolyTs;
                prev_dis = true;
                prev_c2 = c2;
                continue;
            }
        }
        prev_dis = false;
        prev_c2 = c2;

        // create observations for current covered interval
        //std::cout << "   parsed interval: " << c1 << " - " << c2 << " ..."  << std::endl;
        Observations observations;
        observations.truncCounts = infix(contigObservationsF.truncCounts, c1, c2);
        observations.contigId = contigId;
        if (!empty(options.rpkmFileName))
        {
            observations.rpkms = infix(contigCovsF, c1, c2); 
        }
        if (options.useFimoScore)
        {
            observations.fimoScores = infix(contigCovsFimo[0], c1, c2); 
            observations.motifIds = infix(motifIds[0], c1, c2); 
        }
        //std::cout << "TEST: append c1 " << c1 << " c2: " << c2 << std::endl;
        appendValue(data.setObs[0], observations, Generous());
        appendValue(data.setPos[0], c1, Generous());
    }
    // REVERSE 
    unsigned i1_R = length(contigObservationsR.truncCounts) - i2;
    unsigned i2_R = length(contigObservationsR.truncCounts) - i1;
    i = i1_R;
    prev_c1 = i1_R;
    prev_c2 = i1_R;
    prev_dis = false;
    if (options.verbosity >= 2) std::cout << "R: Parse covered intervals and get observations  ..." << "i1_R: " << i1_R << " i2_R: " << i2_R << std::endl;
    while (i < i2_R && (prev_c2 < i2_R && !prev_dis))
    {
        while (i < i2_R && contigObservationsR.truncCounts[i] == 0) ++i;    // find begin of covered interval
        c1 = i;
        if (((int)c1 - (int)options.intervalOffset) > (int)prev_c2)    // if gap bigger than intervalOffset, shift c1 to left
        {
            c1 -= options.intervalOffset;
        }
        else if (prev_c2 == i1_R)           // else if begin of user defined interval. set c1 to this
        {
            c1 = i1_R;
        }
        else if (prev_dis)                  // else if prev. interval was already discarded (and hence not appended)
        {
            c1 = prev_c1;
            ++i;
            if (i >= i2_R) break;
                
            while (i < i2_R && contigObservationsR.truncCounts[i] > 0) ++i;      // find end of covered interval
            c2 = std::min(i + options.intervalOffset, i2_R);
            prev_c2 = c2;
            continue;
        }
        else       
        {                                   // else merge with previous interval 
            c1 = prev_c1;
            eraseBack(data.setObs[1]); 
            eraseBack(data.setPos[1]);
        }
        prev_c1 = c1;
        ++i;

        if (i >= i2_R) break;

        while (i < i2_R && contigObservationsR.truncCounts[i] > 0) ++i;      // find end of covered interval
        c2 = std::min(i + options.intervalOffset, i2_R);

        if (excludePolyA) // check if covered interval contains internal polyA  
        {
            if (checkForPolyT(infix(store.contigStore[contigId].seq, length(store.contigStore[contigId].seq)-c2-1, length(store.contigStore[contigId].seq)-c1-1), options)) 
            {
                ++countPolyAs;
                prev_dis = true;
                prev_c2 = c2;
                continue;
            }
        }
        if (excludePolyT) // check for polyT (polyU)  
        {
            if (checkForPolyA(infix(store.contigStore[contigId].seq, length(store.contigStore[contigId].seq)-c2-1, length(store.contigStore[contigId].seq)-c1-1), options)) 
            {
                ++countPolyTs;
                prev_dis = true;
                prev_c2 = c2;
                continue;
            }
        }
        prev_dis = false;
        prev_c2 = c2;

        // create observations for current covered interval
        Observations observations; 
        observations.truncCounts = infix(contigObservationsR.truncCounts, c1, c2);
        observations.contigId = contigId;
        if (options.useFimoScore)
        {
            observations.fimoScores = infix(contigCovsFimo[1], c1, c2);
            observations.motifIds = infix(motifIds[1], c1, c2); 
        }
        if (!empty(options.rpkmFileName))
        {
            observations.rpkms = infix(contigCovsR, c1, c2);
        }
        appendValue(data.setObs[1], observations, Generous());
        appendValue(data.setPos[1], c1, Generous());
        // corresponding original positions are: t = len - t_R - 1
        // std::cout << "   parsed interval: " << c1 << " - " << c2 << " ..." << std::endl;
    }
    if (options.verbosity >= 2) 
    {
        if (excludePolyA) std::cout << " Excluded " << countPolyAs << " covered intervals from analysis because of internal polyA sites! " << std::endl;
        if (excludePolyT) std::cout << " Excluded " << countPolyTs << " covered intervals from analysis because of internal polyU sites! " << std::endl;
        std::cout << " No. of remaining intervals: " << (length(data.setObs[0]) + length(data.setObs[1])) << "   F: " << length(data.setObs[0]) << "   R: " << length(data.setObs[1]) << std::endl;
    }
    cleanCoveredIntervals(data, length(store.contigStore[contigId].seq), learning, options);
    if (options.verbosity >= 2) 
        std::cout << " No. of remaining intervals after cleaning up: " << (length(data.setObs[0]) + length(data.setObs[1])) << "   F: " << length(data.setObs[0]) << "   R: " << length(data.setObs[1]) << std::endl;
}


// simple linear regression: kde -> N (window count)
template <typename TOptions>
void computeSLR(double &b0, double &b1, Data &data, TOptions &options) // TODO check result
{
    unsigned w_50 = floor((double)options.bandwidthN - 0.1);    // binSize should be odd

    String<double> kdes;
    String<unsigned> counts;
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            // KDE - window truncCount relationship 
            for (unsigned t = 0; t < data.setObs[s][i].length(); ++t)
            {
                unsigned sum = 0;
                for (unsigned j = std::max((int)t - (int)w_50, (int)0); (j < data.setObs[s][i].length()) && (j <= t + w_50); ++j)  // inefficient... update on the fly
                    sum += data.setObs[s][i].truncCounts[j];
                
                appendValue(kdes, data.setObs[s][i].kdesN[t], Generous());
                appendValue(counts, sum, Generous());
                // = std::max(sum, (unsigned)1);  // TODO avoid becoming 0 !
                //out << setObsF[i].kdes[t] << '\t' << sum << '\n';
            }
        }
    }
    std::cout << "  Compute SLR ... " << std::endl;
    double mean_kde = 0.0;
    double mean_count = 0.0;
    for (unsigned i = 0; i < length(kdes); ++i)
    {
        mean_kde += kdes[i];
        mean_count += counts[i];
    }
    mean_kde = mean_kde/(double)length(kdes);
    mean_count = (double)mean_count/(double)length(counts);

    double sum1 = 0;
    double sum2 = 0;
    for (unsigned i = 0; i < length(kdes); ++i)
    {
       sum1 += (kdes[i] - mean_kde) * ((double)counts[i] - mean_count);
       sum2 += pow((kdes[i] - mean_kde), 2);
    }

    b1 = sum1 / sum2;
    b0 = mean_count - b1*mean_kde;

    if (options.verbosity >= 2) std::cout << "Simple linear regression (count <- kde): b0 = " << b0 << " and b1 = " << b1 << " ." << std::endl;
}


template <typename TBai, typename TStore, typename TOptions>
void preproCoveredIntervals(Data &data, double &b0, double &b1, TBai &inputBaiIndex, TStore &store, bool parallelize, TOptions &options)
{
    if (options.verbosity >= 1) std::cout << "  Compute KDEs ... " << std::endl;
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(parallelize ? options.numThreads : 1)) 
#endif
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            data.setObs[s][i].computeKDEs(options);
        }
    } 
    
    if (options.verbosity >= 1) std::cout << "  Estiamte Ns ... " << options.estimateNfromKdes << std::endl;
    if (options.estimateNfromKdes && b0 == 0.0 && b1 == 0.0) 
        computeSLR(b0, b1, data, options);

    // estimate Ns (bin(k; p, N)): either by using raw counts or by using KDEs
    for (unsigned s = 0; s < 2; ++s)
    {
#if HMM_PARALLEL
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(parallelize ? options.numThreads : 1)) 
#endif
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            if (options.estimateNfromKdes) 
                data.setObs[s][i].estimateNs(b0, b1, options);
            else
                data.setObs[s][i].estimateNs(options);

            clear(data.setObs[s][i].kdesN); // only used to estimate Ns
        }
    }

    // if input BAM file given
    if (options.useCov_RPKM && !empty(options.inputBamFileName))
    {
        loadBAMCovariates(data, inputBaiIndex, store, parallelize, options);       // interval-wise
    }
}


template<typename TGAMMA, typename TBIN>
bool learnHMM(Data &data, 
              ModelParams<TGAMMA, TBIN> &modelParams,
              unsigned &contigLen,
              AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    if (options.verbosity >= 1) std::cout << "Build HMM ..." << std::endl;
    HMM<TGAMMA, TBIN> hmm(4, data.setObs, data.setPos, contigLen);       
    hmm.transMatrix = modelParams.transMatrix;
    if (options.verbosity >= 1) 
    {
        myPrint(hmm);
        std::cout << std::endl;
        std::cout << "Baum-Welch  ..." << std::endl;
    }
    CharString learnTag;
    if (options.verbosity >= 1)  std::cout << "            learn binomial parameter" << std::endl;      // in rare cases (with certain parameter combinations) beneficial to learn binomial parameters first
    learnTag = "LEARN_BINOMIAL"; 
    if (!hmm.baumWelch(modelParams, learnTag, options))
        return false;

    if (options.verbosity >= 1)  std::cout << "            learn gamma parameter" << std::endl;
    learnTag = "LEARN_GAMMA";
    if (!hmm.baumWelch(modelParams, learnTag, options))
        return false;

    if (options.verbosity >= 1)  std::cout << "            learn binomial parameter" << std::endl;
    learnTag = "LEARN_BINOMIAL"; 
    if (!hmm.baumWelch(modelParams, learnTag, options))
        return false;

    if (options.verbosity >= 1) std::cout << "            learn gamma parameter" << std::endl;
    learnTag = "LEARN_GAMMA";
    if (!hmm.baumWelch(modelParams, learnTag, options))
        return false;
    // TODO optimize!

    modelParams.transMatrix = hmm.transMatrix;
    data.statePosteriors = hmm.statePosteriors;

    if (options.verbosity >= 2) myPrint(modelParams.gamma1);
    if (options.verbosity >= 2) myPrint(modelParams.gamma2);    
   
#ifdef HMM_PROFILE
    Times::instance().time_learnHMM += (sysTime() - timeStamp);
#endif

    return true;
}


template <typename TGamma, typename TBIN, typename TStore, typename TOptions>
bool learnModel(String<ModelParams<TGamma, TBIN> > &modelParams, 
                String<BamIndex<Bai> > &baiIndices, 
                BamIndex<Bai> &inputBaiIndex, 
                TStore &store, 
                TOptions &options)
{
    resize(baiIndices, length(options.baiFileNames));
    ////////////////////////////////////////////////////
    // for each replicate, learn HMM parameters
    ////////////////////////////////////////////////////
    for (unsigned rep = 0; rep < length(options.baiFileNames); ++rep)
    {
        if (options.verbosity >= 1 && length(options.baiFileNames) == 1) std::cout << "Learn HMM parameters: " << std::endl;
        else if (options.verbosity >= 1 && length(options.baiFileNames) > 1) std::cout << "Learn HMM parameters for replicate: " << rep << std::endl;
        
        // Read BAI index.
        BamIndex<Bai> baiIndex;
        if (!open(baiIndex, toCString(options.baiFileNames[rep])))
        {
            std::cerr << "ERROR: Could not read BAI index file " << options.baiFileNames[rep] << "\n";
            return false;
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
        if (stop) return false;

        // learn KDE - N relationship on all contigs used for other parameter learning as well
        for (unsigned s = 0; s < 2; ++s)
            for (unsigned i = 0; i < length(data.setObs[s]); ++i)
                data.setObs[s][i].computeKDEs(options);

        if (options.estimateNfromKdes) 
            computeSLR(modelParams[rep].slr_NfromKDE_b0, modelParams[rep].slr_NfromKDE_b1, data, options);

        // precompute KDE values, estimate Ns, etc.
        preproCoveredIntervals(data, modelParams[rep].slr_NfromKDE_b0, modelParams[rep].slr_NfromKDE_b1, inputBaiIndex, store, true, options);   

        if (options.verbosity >= 1) std::cout << "Prior ML estimation of density distribution parameters using predefined cutoff ..." << std::endl;
        prior_mle(modelParams[rep].gamma1, modelParams[rep].gamma2, data, options);
        estimateTransitions(modelParams[rep].transMatrix, modelParams[rep].gamma1, modelParams[rep].gamma2, modelParams[rep].bin1, modelParams[rep].bin2, data, options);

        unsigned contigLen = 0; // should not be used within learning
        if (!learnHMM(data, modelParams[rep], contigLen, options))
            return false;

        clear(contigObservationsF);
        clear(contigObservationsR);
        clear(data);
    }
    return true;
}



// in parallel for replicates
template<typename TGAMMA, typename TBIN>
bool applyHMM(Data &newData,
              String<Data> &data_replicates, 
              String<ModelParams<TGAMMA, TBIN> > &modelParams,
              unsigned &contigLen,
              AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    String<HMM<TGAMMA, TBIN> > hmms_replicates;
    reserve(hmms_replicates, length(data_replicates));
    for (unsigned rep = 0; rep < length(data_replicates); ++rep)
    {
        if (options.verbosity >= 1 && length(options.baiFileNames) == 1) std::cout << "Apply learned HMM: " << std::endl;
        else if (options.verbosity >= 1 && length(options.baiFileNames) > 1) std::cout << "Apply learned HMM for replicate: " << rep << std::endl;

        appendValue(hmms_replicates, HMM<TGAMMA, TBIN>(4, data_replicates[rep].setObs, data_replicates[rep].setPos, contigLen));
        auto & hmm = hmms_replicates[rep];
        hmm.transMatrix = modelParams[rep].transMatrix;
        if (options.verbosity >= 1) std::cout << "   applyParameters" << std::endl;
        if (!hmm.applyParameters(modelParams[rep], options))
             return false;

    }

    // merge HMMs and multiply eProbs! (or if not multiple replicates: just using old HMM)
    if (options.verbosity >= 1 && length(options.baiFileNames) > 1) std::cout << "Merge HMMs from replicates ... " << std::endl;
    HMM<TGAMMA, TBIN> mergedHmm = merge_HMMs(hmms_replicates, modelParams, options);

    if (length(options.baiFileNames) > 1)
    {
        // update
        if (options.verbosity >= 1) std::cout << "Update transition probabilities and get final posterior probabilities ... " << std::endl;    
        mergedHmm.updateTransAndPostProbs(options);
    }

    // get states
    mergedHmm.posteriorDecoding(newData.states);

    if (options.useCov_RPKM)    // && !options.g1_k_le_g2_k ?  NOTE: otherwise not necessary, since gamma1.k <= 1
        mergedHmm.rmBoarderArtifacts(newData.states, data_replicates, modelParams);

    newData.statePosteriors = mergedHmm.statePosteriors;
    newData.setObs = mergedHmm.setObs;
    newData.setPos = mergedHmm.setPos;

#ifdef HMM_PROFILE
    Times::instance().time_applyHMM += (sysTime() - timeStamp);
#endif

    return true;
}


template <typename TGamma, typename TBIN, typename TStore, typename TOptions>
bool applyModel(String<String<BedRecord<Bed6> > > &bedRecords_sites, 
                String<String<BedRecord<Bed6> > > &bedRecords_regions, 
                String<ModelParams<TGamma, TBIN> > &modelParams, 
                String<BamIndex<Bai> > &baiIndices, 
                BamIndex<Bai> &inputBaiIndex, 
                TStore &store, 
                TOptions &options)
{
    if (options.verbosity >= 1) std::cout << "Apply learned parameters to whole dataset ..." << std::endl;
    bool stop = false;
    
#if HMM_PARALLEL
    omp_set_num_threads(options.numThreadsA);
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreadsA/length(options.baiFileNames)))    // TODO improve general parallelization concept 
#endif 
    for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
    {
        unsigned contigId = options.applyChr_contigIds[i];
        unsigned contigLen = length(store.contigStore[contigId].seq);
        bool skipContig = false;

        if (options.verbosity >= 1) std::cout << "  " << store.contigNameStore[contigId] << std::endl;

        String<ContigObservations> contigObservationsF;
        String<ContigObservations> contigObservationsR;
        resize(contigObservationsF, length(options.baiFileNames));
        resize(contigObservationsR, length(options.baiFileNames));
        String<Data> data_replicates;
        resize(data_replicates, length(options.baiFileNames));    

        for (unsigned rep = 0; rep < length(options.baiFileNames); ++rep)
        {
            if (options.verbosity >= 1 && length(options.baiFileNames) == 1) std::cout << "Get preprocessed intervals: " << std::endl;
            else if (options.verbosity >= 1 && length(options.baiFileNames) > 1) std::cout << "Get preprocessed intervals for replicate: " << rep << std::endl;
            

            // for each replicate get preprocessed covered intervals
            int r = loadObservations(contigObservationsF[rep], contigObservationsR[rep], contigId, options.bamFileNames[rep], baiIndices[rep], store, options);
            if (r == 1) // error
            {
                SEQAN_OMP_PRAGMA(critical)
                stop = true; 
            }
            else if (r == 2)    // no alignments (F & R) for one replicate, ignore contig
            {
                SEQAN_OMP_PRAGMA(critical)
                skipContig = true; 
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
                    preproCoveredIntervals(c_data, modelParams[rep].slr_NfromKDE_b0, modelParams[rep].slr_NfromKDE_b1, inputBaiIndex, store, false, options);
                    data_replicates[rep] = c_data;
                }
                else
                {
                    // if no covered regions remaining for one replicate, ignore contig
                    SEQAN_OMP_PRAGMA(critical)
                    skipContig = true;
                }
            }
        }
        if (stop || skipContig) continue;

        // get intersection of covered intervals and clip individual observations
        if (length(options.baiFileNames) > 1)
        {
            if (options.verbosity >= 1) std::cout << "Intersect covered intervals from replicates ... (contig: " << store.contigNameStore[contigId] << ")." << std::endl;
            intersect_replicateIntervals(data_replicates);
        }
        
        Data newData = data_replicates[0];
        // build individual HMMs on clipped intervals, merge HMMs, update trans. probs, compute post. probs and get states
        if (!applyHMM(newData, data_replicates, modelParams, contigLen, options))
        {
            stop = true;
            continue;
        }

        writeStates(bedRecords_sites[i], newData, store, contigId, options);  
        if (!empty(options.outRegionsFileName))
        {
            writeRegions(bedRecords_regions[i], newData, store, contigId, options);              
        }
    }
    if (stop) return false;

    return true;
}


template <typename TGamma, typename TBIN, typename TOptions>
bool doIt(String<ModelParams<TGamma, TBIN> > &modelParams, TOptions &options)
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

    setSomeParameters(modelParams, options);

    //////////////////////////////////////////////////////////////////

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

    // learn model
    String<BamIndex<Bai> > baiIndices;    
    if (!learnModel(modelParams, baiIndices, inputBaiIndex, store, options))
        return 1;


    // store contig-wise bedRecords
    String<String<BedRecord<Bed6> > > bedRecords_sites;
    String<String<BedRecord<Bed6> > > bedRecords_regions;
    resize(bedRecords_sites, length(options.applyChr_contigIds), Exact());
    resize(bedRecords_regions, length(options.applyChr_contigIds), Exact());
    // apply model to whole dataset
    if (!applyModel(bedRecords_sites, bedRecords_regions, modelParams, baiIndices, inputBaiIndex, store, options))
        return 1;

    if (options.verbosity >= 2) std::cout << "Write bedRecords to BED file ... " << options.outFileName << std::endl;
    // write crosslink sites
    BedFileOut outBed(toCString(options.outFileName)); 
    for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
    {
        for (unsigned j = 0; j < length(bedRecords_sites[i]); ++j)
            writeRecord(outBed, bedRecords_sites[i][j]);
    }
    // write binding regions
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
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
    //std::cout << "  Time needed for loadObservations: " << Times::instance().time_loadObservations/60.0 << "min" << std::endl;
    //std::cout << "  Time needed for learnHMM: " << Times::instance().time_learnHMM/60.0 << "min" << std::endl;
    //std::cout << "  Time needed for applyHMM: " << Times::instance().time_applyHMM/60.0 << "min" << std::endl;
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
        if (length(options.baiFileNames) == 1) out << "Parameter learned: " << std::endl;
        else out << "Parameter learned for replicate: " << rep << std::endl;
        printParams(out, modelParams[rep].gamma1, 1);
        printParams(out, modelParams[rep].gamma2, 2);
        printParams(out, modelParams[rep].bin1, 1);
        printParams(out, modelParams[rep].bin2, 2);
        printParams(out, modelParams[rep].transMatrix);
        out << std::endl;
        out << "slr_NfromKDE.b0" << '\t' << modelParams[rep].slr_NfromKDE_b0 << std::endl;
        out << "slr_NfromKDE.b1" << '\t' << modelParams[rep].slr_NfromKDE_b1 << std::endl;  
    }
    return 0;
}

    
#endif
