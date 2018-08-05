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



template <typename TContigObservations, typename TStore>
int loadObservations(TContigObservations &contigObservationsF, TContigObservations &contigObservationsR, unsigned contigId, TStore &store, AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    if (options.verbosity >= 2) std::cout << "Parse alignments ... " << std::endl;

    // Open BamFileIn for reading.
    if (options.verbosity >= 2) std::cout << "Open Bam and Bai file ... "  << "\n";
    BamFileIn inFile;
    if (!open(inFile, toCString(options.bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << options.bamFileName << " for reading.\n";
        return 1;
    }
    BamHeader header;
    readHeader(header, inFile);
    // Read BAI index.
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(options.baiFileName)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << options.baiFileName << "\n";
        return 1;
    }
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


template <typename TStore, typename TOptions>
bool loadBAMCovariates(Data &data, TStore &store, TOptions &options)
{
    if (options.verbosity >= 2) std::cout << "Parse alignments ... " << std::endl;

    // Open BamFileIn for reading.
    std::cout << "Open Bam and Bai file ... "  << "\n";
    BamFileIn inFile;
    if (!open(inFile, toCString(options.inputBamFileName)))
    {
        std::cerr << "ERROR: Could not open " << options.inputBamFileName << " for reading.\n";
        return false;
    }
    BamHeader header;
    readHeader(header, inFile);
    // Read BAI index.
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(options.inputBaiFileName)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << options.inputBaiFileName << "\n";
        return false;
    }

    if (options.verbosity >= 1) std::cout << "  Parse input BAM, get truncCounts, compute KDEs ... " << std::endl;
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            // Translate from contig name to rID.
            int rID = 0;
            if (!getIdByName(rID, contigNamesCache(context(inFile)), store.contigNameStore[data.setObs[s][i].contigId]))
            {
                std::cerr << "ERROR: Contig " << store.contigNameStore[data.setObs[s][i].contigId] << " not known.\n";
                return false; 
            }
            String<__uint16> truncCounts;
            resize(truncCounts, data.setObs[s][i].length(), 0, Exact());
            if (s == 0)
            {
                unsigned beginPos = data.setPos[s][i];
                unsigned endPos = data.setPos[s][i] + data.setObs[s][i].length();
                parse_bamRegion(truncCounts, inFile, baiIndex, rID, beginPos, endPos, true, options);
            }
            else
            {
                unsigned beginPos = length(store.contigStore[data.setObs[s][i].contigId].seq) - (data.setObs[s][i].length() + data.setPos[s][i]);
                unsigned endPos = beginPos + data.setObs[s][i].length();

                parse_bamRegion(truncCounts, inFile, baiIndex, rID, beginPos, endPos, false, options);
                // reverse
                reverse(truncCounts);
            }
            // compute KDEs
            data.setObs[s][i].computeKDEs(truncCounts, options);
        }
    }

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
    while (i < i2)
    {

        while (i < i2 && contigObservationsF.truncCounts[i] == 0) ++i;    // find begin of covered interval     
        c1 = i;
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
    while (i < i2_R)
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
        std::cout << " Excluded " << countPolyAs << " covered intervals from analysis because of internal polyA sites! " << std::endl;
        std::cout << " Excluded " << countPolyTs << " covered intervals from analysis because of internal polyU sites! " << std::endl;
        std::cout << " No. of remaining intervals: " << (length(data.setObs[0]) + length(data.setObs[1])) << "   F: " << length(data.setObs[0]) << "   R: " << length(data.setObs[1]) << std::endl;
    }
    cleanCoveredIntervals(data, length(store.contigStore[contigId].seq), learning, options);
    if (options.verbosity >= 2) 
        std::cout << " No. of remaining intervals after cleaning up: " << (length(data.setObs[0]) + length(data.setObs[1])) << std::endl;
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


template <typename TStore, typename TOptions>
void preproCoveredIntervals(Data &data, double &b0, double &b1, TStore &store, TOptions &options)
{
    if (options.verbosity >= 1) std::cout << "  Compute KDEs ... " << std::endl;
    for (unsigned s = 0; s < 2; ++s)
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
            data.setObs[s][i].computeKDEs(options);
    
    
    if (options.verbosity >= 1) std::cout << "  Estiamte Ns ... " << options.estimateNfromKdes << std::endl;
    if (options.estimateNfromKdes && b0 == 0.0 && b1 == 0.0) 
        computeSLR(b0, b1, data, options);

    // estimate Ns (bin(k; p, N)): either by using raw counts or by using KDEs
    for (unsigned s = 0; s < 2; ++s)
    {
        for (unsigned i = 0; i < length(data.setObs[s]); ++i)
        {
            if (options.estimateNfromKdes) 
                data.setObs[s][i].estimateNs(b0, b1, options);
            else
                data.setObs[s][i].estimateNs(options);

            clear(data.setObs[s][i].kdesN);
        }
    }

    // if input BAM file given
    if (options.useCov_RPKM && !empty(options.inputBamFileName))
    {
        loadBAMCovariates(data, store, options);       // interval-wise
    }
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool learnHMM(Data &data, 
              String<String<double> > &transMatrix_1,
              TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2,
              TDOUBLE /**/,
              unsigned &contigLen,
              AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    if (options.verbosity >= 1) std::cout << "Build HMM ..." << std::endl;
    HMM<TGAMMA, TBIN, TDOUBLE> hmm(4, data.setObs, data.setPos, contigLen);       
    hmm.transMatrix = transMatrix_1;
    if (options.verbosity >= 1) 
    {
        myPrint(hmm);
        std::cout << std::endl;
        std::cout << "Baum-Welch  ..." << std::endl;
    }
    if (options.verbosity >= 1)  std::cout << "            learn binomial parameter" << std::endl;
    CharString learnTag = "LEARN_BINOMIAL"; 
    if (!hmm.baumWelch(d1, d2, bin1, bin2, learnTag, options))
        return false;

    if (options.verbosity >= 1)  std::cout << "            learn gamma parameter" << std::endl;
    learnTag = "LEARN_GAMMA";
    if (!hmm.baumWelch(d1, d2, bin1, bin2, learnTag, options))
        return false;

    if (options.verbosity >= 1)  std::cout << "            learn binomial parameter" << std::endl;
    learnTag = "LEARN_BINOMIAL"; 
    if (!hmm.baumWelch(d1, d2, bin1, bin2, learnTag, options))
        return false;

    if (options.verbosity >= 1) std::cout << "            learn gamma parameter" << std::endl;
    learnTag = "LEARN_GAMMA";
    if (!hmm.baumWelch(d1, d2, bin1, bin2, learnTag, options))
        return false;


    transMatrix_1 = hmm.transMatrix;

    if (options.verbosity >= 2) myPrint(d1);
    if (options.verbosity >= 2) myPrint(d2);
    if (options.posteriorDecoding)
        hmm.posteriorDecoding(data.states);
    else
        hmm.viterbi_log(data.states);
    
    if (options.useCov_RPKM && !options.g1_k_le_g2_k)    // NOTE: otherwise not necessary, since gamma1.k <= 1 or gamma1.k <= gamma2.k
        hmm.rmBoarderArtifacts(data.states, d1);
    data.statePosteriors = hmm.statePosteriors;
   
#ifdef HMM_PROFILE
    Times::instance().time_learnHMM += (sysTime() - timeStamp);
#endif

    return true;
}


template<typename TGAMMA, typename TBIN, typename TDOUBLE>
bool applyHMM(Data &data, 
              String<String<double> > &transMatrix_1,
              TGAMMA &d1, TGAMMA &d2, TBIN &bin1, TBIN &bin2,
              TDOUBLE /**/,
              unsigned &contigLen,
              AppOptions &options)
{
#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif
    if (options.verbosity >= 1) std::cout << "   build HMM" << std::endl;
    HMM<TGAMMA, TBIN, TDOUBLE> hmm(4, data.setObs, data.setPos, contigLen);
    hmm.transMatrix = transMatrix_1;
    if (!hmm.applyParameters(d1, d2, bin1, bin2, options))
        return false;

    if (options.posteriorDecoding)
        hmm.posteriorDecoding(data.states);
    else
        hmm.viterbi_log(data.states);

    if (options.verbosity >= 2)
    {
        std::cout << "Intermediate parameters after applyHMM():" << std::endl;
        myPrint(hmm);
        myPrint(d1);
        myPrint(d2);
    }
 
    if (options.useCov_RPKM)    // NOTE: otherwise not necessary, since gamma1.k <= 1
        hmm.rmBoarderArtifacts(data.states, d1);
    data.statePosteriors = hmm.statePosteriors;

#ifdef HMM_PROFILE
    Times::instance().time_applyHMM += (sysTime() - timeStamp);
#endif

    return true;
}


inline
const char * myTempFileName(std::string const suffix = "", std::string const tempPath = "")
{
    if (mkdir(tempPath.c_str(), 0777) == -1)
    {
        if(errno == EEXIST ) {
            //std::cout << "Directory already exists " << std::endl;
        } else {
            std::cerr << tempPath << std::endl;
            throw std::runtime_error("ERROR: Could not create directory");
        }
    }

    static char fileNameBuffer[1000];
    strcpy(fileNameBuffer, (tempPath + "PURECLIP.XXXXXX" + suffix).c_str());

    int _tmp = mkstemps(fileNameBuffer, suffix.size());
    if (_tmp == -1)
        throw std::runtime_error("Could not create temp file");

    close(_tmp);    
    return fileNameBuffer;
}

inline bool exists_test(const CharString& fileName) {
    std::ifstream f(toCString(fileName));
    return f.good();
}


template <typename TGamma, typename TBIN, typename TDOUBLE, typename TOptions>
bool doIt(TGamma &gamma1, TGamma &gamma2, TBIN &bin1, TBIN &bin2, TDOUBLE /**/, TOptions &options)
{
#if HMM_PARALLEL
    omp_set_num_threads(options.numThreads);
#endif  

#ifdef HMM_PROFILE
    double timeStamp = sysTime();
#endif

    typedef  FragmentStore<>    TStore;
    TStore store;
    if (options.verbosity >= 1) std::cout << "Loading reference ... " << std::endl;;
    
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
    if (!options.useHighPrecision) 
    {
        options.min_eProbSum = pow((double)10.0, DBL_MIN_10_EXP + 50);
        std::cout << " Set options.min_eProbSum : " << options.min_eProbSum << std::endl;  
    }
    else 
    {
        options.min_eProbSum = pow((long double)10.0, LDBL_MIN_10_EXP + 100);
        std::cout << " Set options.min_eProbSum : " << options.min_eProbSum << std::endl;
    }    
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

    // *****************
    double slr_NfromKDE_b0 = 0.0;
    double slr_NfromKDE_b1 = 0.0;  


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

        int r = loadObservations(contigObservationsF[i], contigObservationsR[i], contigId, store, options);
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
        computeSLR(slr_NfromKDE_b0, slr_NfromKDE_b1, data, options);
  
    // precompute KDE values, estimate Ns, etc.
    preproCoveredIntervals(data, slr_NfromKDE_b0, slr_NfromKDE_b1, store, options);
 
    gamma1.tp = options.useKdeThreshold;
    gamma2.tp = options.useKdeThreshold;       // left tuncated    

    if (options.verbosity >= 1) std::cout << "Prior ML estimation of density distribution parameters using predefined cutoff ..." << std::endl;
    prior_mle(gamma1, gamma2, data, options);
    String<String<double> > transMatrix;
    estimateTransitions(transMatrix, gamma1, gamma2, bin1, bin2, data, options);

    unsigned contigLen = 0; // should not be used within learning
    if (!learnHMM(data, transMatrix, gamma1, gamma2, bin1, bin2, (TDOUBLE)0.0, contigLen, options))
        return 1;

    clear(contigObservationsF);
    clear(contigObservationsR);

    if (options.verbosity >= 1) std::cout << "Apply learned parameters to whole genome  ..." << std::endl;
#if HMM_PARALLEL
    omp_set_num_threads(options.numThreadsA);
#endif  
    String<CharString> contigTempFileNamesBed;
    String<CharString> contigTempFileNamesBed2;  
    resize(contigTempFileNamesBed, length(store.contigStore));
    resize(contigTempFileNamesBed2, length(store.contigStore));

#ifdef HMM_PROFILE
    double timeStamp2 = sysTime();
#endif

#if HMM_PARALLEL
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreadsA)) 
#endif  
    for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
    {
        unsigned contigId = options.applyChr_contigIds[i];
        if (options.verbosity >= 1) std::cout << "  " << store.contigNameStore[contigId] << std::endl;

        ContigObservations contigObservationsF;
        ContigObservations contigObservationsR;

        int r = loadObservations(contigObservationsF, contigObservationsR, contigId, store, options);
        if (r == 1)
        {
            stop = true; 
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
            extractCoveredIntervals(c_data, contigObservationsF, contigObservationsR, c_contigCovsF, c_contigCovsR, c_contigCovsFimo, c_motifIds, contigId, i1, i2, options.excludePolyA, options.excludePolyT, false, store, options); 

            if (!empty(c_data.setObs[0]) || !empty(c_data.setObs[1]))   // TODO handle cases
            {

                preproCoveredIntervals(c_data, slr_NfromKDE_b0, slr_NfromKDE_b1, store, options);

                // Apply learned parameters
                unsigned contigLen = length(store.contigStore[contigId].seq);
                if (!applyHMM(c_data, transMatrix, gamma1, gamma2, bin1, bin2, (TDOUBLE)0.0, contigLen, options))
                {
                    SEQAN_OMP_PRAGMA(critical)
                    stop = true;
                }

                // Temp. output
                CharString tempFileNameBed;

                if (empty(options.tempPath))
                {
                    SEQAN_OMP_PRAGMA(critical)
                    append(tempFileNameBed, SEQAN_TEMP_FILENAME());
                    std::stringstream ss;
                    ss << contigId;
                    append(tempFileNameBed, ss.str());
                    append(tempFileNameBed, ".bed");
                }
                else
                {
                    SEQAN_OMP_PRAGMA(critical)
                    tempFileNameBed = myTempFileName(".bed", toCString(options.tempPath));
                }

                contigTempFileNamesBed[contigId] = tempFileNameBed;
                if (options.verbosity >= 2) std::cout << "temp file Name: " << tempFileNameBed << std::endl;
                BedFileOut outBed(toCString(tempFileNameBed)); 
                writeStates(outBed, c_data, store, contigId, options);  
                
                if (!empty(options.outRegionsFileName))
                {
                    if (empty(options.tempPath))
                    {

                        SEQAN_OMP_PRAGMA(critical)
                        tempFileNameBed = SEQAN_TEMP_FILENAME();
                        std::stringstream ss;
                        ss << contigId;
                        append(tempFileNameBed, ss.str());
                        append(tempFileNameBed, ".regions.bed");
                    }
                    else
                    {
                        SEQAN_OMP_PRAGMA(critical)
                        tempFileNameBed = myTempFileName(".regions.bed", toCString(options.tempPath));
                    }

                    contigTempFileNamesBed2[contigId] = tempFileNameBed;
                    if (options.verbosity >= 2) std::cout << "temp file Name: " << tempFileNameBed << std::endl;
                    BedFileOut outBed2(toCString(tempFileNameBed)); 

                    writeRegions(outBed2, c_data, store, contigId, options);              
                }
            }
        }
    }
    if (stop) return 1;

    // Append content of temp files to final output in contig order
    // crosslink sites
    BedFileOut outBed(toCString(options.outFileName)); 
    for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
    {
        unsigned contigId = options.applyChr_contigIds[i];
        if (!empty(contigTempFileNamesBed[contigId]))
        {
            // BED
            BedFileIn bedFileIn;
            if (!open(bedFileIn, toCString(contigTempFileNamesBed[contigId])))
            {
                //std::cerr << "ERROR: Could not open temporary bed file: " << contigTempFileNamesBed[contigId] << "\n";
                continue;
            }

            BedRecord<seqan::Bed6> bedRecord;
            while (!atEnd(bedFileIn))
            {
                readRecord(bedRecord, bedFileIn);
                writeRecord(outBed, bedRecord);
            }
            std::remove(toCString(contigTempFileNamesBed[contigId]));
            if (exists_test(contigTempFileNamesBed[contigId]))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << contigTempFileNamesBed[contigId]  << std::endl;
                return 1;
            }
        }
    }

    // Append content of temp files to final output in contig order
    // binding regions
    if (!empty(options.outRegionsFileName))
    {
        BedFileOut outBed2(toCString(options.outRegionsFileName)); 
        for (unsigned i = 0; i < length(options.applyChr_contigIds); ++i)
        {
            unsigned contigId = options.applyChr_contigIds[i];
            if (!empty(contigTempFileNamesBed2[contigId]))
            {
                // BED
                BedFileIn bedFileIn;
                if (!open(bedFileIn, toCString(contigTempFileNamesBed2[contigId])))
                {
                    //std::cerr << "ERROR: Could not open temporary bed file: " << contigTempFileNamesBed[contigId] << "\n";
                    continue;
                }

                BedRecord<seqan::Bed6> bedRecord;
                while (!atEnd(bedFileIn))
                {
                    readRecord(bedRecord, bedFileIn);
                    writeRecord(outBed2, bedRecord);
                }
                std::remove(toCString(contigTempFileNamesBed2[contigId]));
                if (exists_test(contigTempFileNamesBed2[contigId]))
                {
                    std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << contigTempFileNamesBed2[contigId]  << std::endl;
                    return 1;
                }
            }
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

    CharString fileNameStats = options.outFileName;
    append(fileNameStats, ".stats");
    std::ofstream out(toCString(fileNameStats), std::ios::binary | std::ios::out);
    printParams(out, gamma1, 1);
    printParams(out, gamma2, 2);
    out << "options.useKdeThreshold" << '\t' << options.useKdeThreshold << std::endl;

    return 0;
}

 

// transcript wise?
    
#endif
