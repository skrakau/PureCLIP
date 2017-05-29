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

#include <fstream>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h> 

#include <seqan/arg_parse.h>

#include <iostream>
#include <seqan/find.h>

#include <seqan/bed_io.h>
#include <seqan/misc/interval_tree.h> 

using namespace seqan;


struct AppOptions
{
    CharString inRefFileName;
    CharString inBamFileName;
    CharString inBedFileName;
    CharString outFileName;
    CharString sitesFileName;

    unsigned windowSize;
    bool useGivenWindow;
    bool addScoreToName;

    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;


    AppOptions() :
        windowSize(100),
        useGivenWindow(false),
        addScoreToName(false),  
        verbosity(1)
    {}
};


ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("winextract");
    // Set short description, version, and date.
    setShortDescription(parser, "Extract sequence window");
    setVersion(parser, "0.1");
    setDate(parser, "Mai 2017");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "Extracts sequences form reference given BED regions.");

    // We require one argument.
    addOption(parser, ArgParseOption("g", "genome", "Genome reference file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "genome", ".fa .fasta");
    setRequired(parser, "genome", true);   
    addOption(parser, ArgParseOption("c", "in-bed", "Input cand-regions.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "in-bed", ".bed");
    setRequired(parser, "in-bed", true);

    addOption(parser, ArgParseOption("o", "output", "Output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output", ".fasta .fa");
    setRequired(parser, "output", true);
    addOption(parser, ArgParseOption("w", "window", "Window size to analyse.", ArgParseArgument::INTEGER));
    setMinValue(parser, "window", "5");
    setMaxValue(parser, "window", "2000"); // ?
     
    addOption(parser, ArgParseOption("u", "uow", "Use given window."));
    addOption(parser, ArgParseOption("s", "asn", "Add score to output sequence name."));

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.inRefFileName, parser, "genome");
    getOptionValue(options.inBedFileName, parser, "in-bed");
    getOptionValue(options.outFileName, parser, "output");
    getOptionValue(options.windowSize, parser, "window");
    
    if (isSet(parser, "uow"))
        options.useGivenWindow = true;

    if (isSet(parser, "asn"))
        options.addScoreToName = true;

    return ArgumentParser::PARSE_OK;
}


template<typename TOptions>
bool doIt(TOptions & options)
{ 
    typedef FragmentStore<> TFragmentStore;

    std::cout << "Load reference file ..." << std::endl;
    TFragmentStore store;
    loadContigs(store, toCString(options.inRefFileName));
   
    std::cout << "Load candidate sites from BED file ..." << std::endl;
    String<String<IntervalAndCargo<> > > candRegs;
    resize(candRegs, length(store.contigNameStore));
    BedFileIn bedIn(toCString(options.inBedFileName));
    BedRecord<Bed6> bedRecord;


    std::cout << "Extract sequence for windows around candidate sites ..." << std::endl;
    SeqFileOut seqFileOut(toCString(options.outFileName));
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
        for (unsigned i = 0; i < length(store.contigNameStore); ++i)
        {
            if (bedRecord.ref == store.contigNameStore[i]) // || bedRecord.ref == prefix(store.contigNameStore[i], length(bedRecord.ref)) )
            {
                unsigned beginPos;
                unsigned endPos;

                if (!options.useGivenWindow)
                {
                    int pos = bedRecord.beginPos;
                    beginPos = 0;
                    endPos = length(store.contigStore[i].seq);

                    if (pos - (int)(options.windowSize/2) > 0) 
                        beginPos = pos - (options.windowSize/2);
                    else 
                        continue;                                                   // for simplicity ignore windows overlapping with contig border
                    if (pos + (options.windowSize/2) + 1 < length(store.contigStore[i].seq) )
                        endPos = pos + (options.windowSize/2) +1;
                    else 
                        continue;
                }
                else
                {
                    beginPos = std::max(bedRecord.beginPos, 0);
                    endPos = std::min(bedRecord.endPos, (int)length(store.contigStore[i].seq));
                }
        
                //std::cout << store.contigNameStore[i] << "  pos: " << pos << "  beginPos: " << beginPos << "  endPos: " << endPos << "  length contig: " << length(store.contigStore[i].seq) << std::endl;
                DnaString inf = infix(store.contigStore[i].seq, beginPos, endPos);
                if (bedRecord.strand == '-') reverseComplement(inf);
                
                CharString id = "cand_";
                std::stringstream ss;
                ss <<  store.contigNameStore[i];
                ss << "_";
                ss << beginPos;
                ss << "_";
                ss << endPos;
                ss << "_";
                if (bedRecord.strand == '+')
                    ss << "F";
                else
                    ss << "R";
                if (options.addScoreToName)
                    {
                    ss << "_";
                    ss << bedRecord.score;
                }

                CharString str = ss.str();  
                ss.str("");  
                ss.clear();  
                append(id, str);

                writeRecord(seqFileOut, id, inf);

                continue;
            }
        }
    }
   
	return 0;
}

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    doIt(options);

    return 0;
}
