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

#ifndef APPS_HMMS_PARSE_ALIGNMENTS_H_
#define APPS_HMMS_PARSE_ALIGNMENTS_H_

#include <iostream>
#include <fstream>

using namespace seqan;



template <typename TContigObservations, typename TBamIn, typename TBai, typename TOptions>
bool parse_bamRegion(TContigObservations &contigObservationsF, TContigObservations &contigObservationsR,  TBamIn &inFile, TBai &baiIndex, int const& rID, TOptions &options)
{
    if (options.verbosity >= 2)
        std::cout << "Parse BAM region " << std::endl;
 
    int jump_beginPos = 0;
    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inFile, hasAlignments, rID, jump_beginPos, length(contigObservationsF.truncCounts)-1, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << jump_beginPos << ":" << (length(contigObservationsF.truncCounts)-1) << "\n";
        return false;
    }
    if (!hasAlignments)
    {
        std::cout << "WARNING: no alignments here " << jump_beginPos << ":" << (length(contigObservationsF.truncCounts)-1) << "\n";
        return false;  
    }

    // Seek linearly to the selected position
    BamAlignmentRecord bamRecord;
    BamAlignmentRecord newBamRecord;
    while (!atEnd(inFile))
    {
        readRecord(bamRecord, inFile);

        // If we are on the next reference 
        if (bamRecord.rID == -1 || bamRecord.rID > rID)
            break;

        if (!hasFlagRC(bamRecord))          // Forward
        {
            if (contigObservationsF.truncCounts[bamRecord.beginPos] < 254)      // uint8, discard interval if > anyway ...
                ++contigObservationsF.truncCounts[bamRecord.beginPos]; 
        }
        else                                // Reverse  
        {
            if (contigObservationsR.truncCounts[bamRecord.beginPos + getAlignmentLengthInRef(bamRecord) - 1] < 254)
                ++contigObservationsR.truncCounts[bamRecord.beginPos + getAlignmentLengthInRef(bamRecord) - 1]; 
        }
    }
    return true;
}




template <typename TTruncCounts, typename TBamIn, typename TBai, typename TOptions>
bool parse_bamRegion(TTruncCounts &truncCounts, TBamIn &inFile, TBai &baiIndex, int const& rID, unsigned beginPos, unsigned endPos, bool isForward, TOptions &options)
{
    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    int jump_beginPos;
    if (isForward)
        jump_beginPos = beginPos;
    else
        jump_beginPos = beginPos - 100;

    if (!jumpToRegion(inFile, hasAlignments, rID, jump_beginPos, (endPos+1000), baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
        return false;
    }
    /*if (!hasAlignments)
    {
        std::cout << "WARNING: no input alignments here " << beginPos << ":" << endPos << " , " << rID << "\n";
        return false;  
    }*/

    // Seek linearly to the selected position
    BamAlignmentRecord bamRecord;
    BamAlignmentRecord newBamRecord;
    while (!atEnd(inFile))
    {
        readRecord(bamRecord, inFile);

        // If we are on the next reference 
        if (bamRecord.rID == -1 || bamRecord.rID > rID)
            break;
        if (bamRecord.beginPos >= (endPos + 100))
            break;

        if (isForward && !hasFlagRC(bamRecord))          // Forward
        {
            if (bamRecord.beginPos >= beginPos &&                       // already there ?
                bamRecord.beginPos < endPos &&                          // before end of interval ?        
                truncCounts[bamRecord.beginPos - beginPos] < 254)      // uint8,  ?
            {
                ++truncCounts[bamRecord.beginPos - beginPos]; 
            }
        }
        else if (!isForward && hasFlagRC(bamRecord))                             // Reverse  
        {
            if ((bamRecord.beginPos + getAlignmentLengthInRef(bamRecord) - 1) >= beginPos &&                        // already there ?
                (bamRecord.beginPos + getAlignmentLengthInRef(bamRecord) - 1) < endPos &&                          // before end of interval  ? 
                 truncCounts[bamRecord.beginPos - beginPos + getAlignmentLengthInRef(bamRecord) - 1] < 254)
            {
                ++truncCounts[bamRecord.beginPos - beginPos + getAlignmentLengthInRef(bamRecord) - 1]; 
            }
        }
    }

    return true;
}


    
#endif
