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
// ==========================================================================
// Author: Sabrina Krakau <krakau@molgen.mpg.de>
// ==========================================================================

#define HMM_PROFILE


#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include <seqan/misc/name_store_cache.h>
#include <seqan/arg_parse.h>

#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>


#include "util.h"
#include "call_sites.h"

using namespace seqan;

 

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("pureclip");
    // Set short description, version, and date.
    setShortDescription(parser, "Protein-RNA crosslink site detection ");
    setVersion(parser, "0.1");
    setDate(parser, "Mai 2017");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <-i \\fIBAM FILE\\fP> <-bai \\fIBAI FILE\\fP> <-g \\fIGENOME FILE\\fP> <-o \\fIOUTPUT BED FILE\\fP> ");
    addDescription(parser, "Protein-RNA crosslink site detection using a non-homogeneous HMM.");

    // We require one argument.
    addOption(parser, ArgParseOption("i", "in", "Target bam file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "in", ".bam");
    setRequired(parser, "in", true);

    addOption(parser, ArgParseOption("bai", "bai", "Target bam index file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "bai", ".bai");
    setRequired(parser, "bai", true);

    addOption(parser, ArgParseOption("g", "genome", "Genome reference file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "genome", ".fa .fasta");
    setRequired(parser, "genome", true);  

    addOption(parser, ArgParseOption("o", "out", "Output file to write crosslink sites.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "out", ".bed");
    setRequired(parser, "out", true);
    addOption(parser, ArgParseOption("or", "or", "Output file to write binding regions.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "or", ".bed");
    addOption(parser, ArgParseOption("p", "par", "Output file to write learned parameters.", ArgParseArgument::OUTPUT_FILE));
    //setRequired(parser, "par", true);
    

    addSection(parser, "Options");

    addOption(parser, ArgParseOption("iv", "inter", "Genomic chromosomes to learn HMM parameters, e.g. 'chr1;chr2;chr3'. Contigs have to be in the same order as in BAM file. Useful to reduce runtime and memory consumption. Default: all contigs from reference file are used (useful when applying to transcript-wise alignments or poor data).", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("chr", "chr", "Contigs to apply HMM, e.g. 'chr1;chr2;chr3;'. Contigs have to be in the same order as in BAM file.", ArgParseArgument::STRING));

    addOption(parser, ArgParseOption("bw", "bdw", "Bandwidth for kernel density estimation. NOTE: Increasing the bandwidth increases runtime and memory consumption. Default: 50.", ArgParseArgument::INTEGER));
    setMinValue(parser, "bdw", "1");
    setMaxValue(parser, "bdw", "500"); 

    addOption(parser, ArgParseOption("dm", "dm", "Distance used to merge individual crosslink sites to binding regions. Default: 8", ArgParseArgument::INTEGER));


    addSection(parser, "Options for incorporating covariates");

    addOption(parser, ArgParseOption("is", "is", "Covariates file: position-wise values, e.g. smoothed reads start counts (KDEs) from input data. ", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "is", ".bed");
    addOption(parser, ArgParseOption("ibam", "ibam", "File containing mapped reads from control experiment, e.g. eCLIP input.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ibam", ".bam");
    addOption(parser, ArgParseOption("ibai", "ibai", "File containing BAM index corresponding to mapped reads from control experiment", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ibai", ".bai");

    addOption(parser, ArgParseOption("fis", "fis", "Fimo input motif score covariates file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "fis", ".bed");
    addOption(parser, ArgParseOption("nim", "nim", "Max. motif ID to use. Default: Only covariates with motif ID 1 are used.", ArgParseArgument::INTEGER));



    addSection(parser, "Advanced user options");

    addOption(parser, ArgParseOption("vtb", "vtb", "Use Viterbi instead of posterior decoding."));
    addOption(parser, ArgParseOption("m", "mibr", "Maximum number of iterations within BRENT algorithm.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mibr", "1");
    setMaxValue(parser, "mibr", "1000");
    addOption(parser, ArgParseOption("w", "mibw", "Maximum number of iterations within Baum-Welch algorithm.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mibw", "0");
    setMaxValue(parser, "mibw", "500");
    addOption(parser, ArgParseOption("g1kmin", "g1kmin", "Minimum shape k of 'non-enriched' gamma distribution (g1.k).", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g1kmax", "g1kmax", "Maximum shape k of 'non-enriched' gamma distribution (g1.k).", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g2kmin", "g2kmin", "Minimum shape k of 'enriched' gamma distribution (g2.k).", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("g2kmax", "g2kmax", "Maximum shape k of 'enriched' gamma distribution (g2.k).", ArgParseArgument::DOUBLE));
    //addOption(parser, ArgParseOption("g1g2k", "g1g2k", "Force 'non-enriched' gamma parameter k <= 'enriched' gamma parameter k."));

    addOption(parser, ArgParseOption("mkn", "mkn", "Max. k/N ratio (read start sites/N) used to learn truncation probabilities for 'non-crosslink' and 'crosslink' emission probabilities (high ratios might originate from mapping artifacts that can disturb parameter learning). Default: 1.0", ArgParseArgument::DOUBLE));
    setMinValue(parser, "mkn", "0.5");
    setMaxValue(parser, "mkn", "1.5"); 

    addOption(parser, ArgParseOption("mtp", "mtp", "Min. transition probability from state '2' to '3' (for poor data, where no clear distinction between 'enriched' and 'non-enriched' is possible). Default: 0.0001.", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("mk", "mkde", "Minimum KDE value used for fitting left-truncated gamma distributions. Default: corresponding to singleton read start.", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("ntp", "ntp", "Only sites with n >= ntp are used to learn binomial probability parameters (bin1.p, bin2.p). Default: 10", ArgParseArgument::DOUBLE));

    addOption(parser, ArgParseOption("pa", "pat", "Length threshold for internal poly-X stretches to get excluded.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ea1", "epal", "Exclude intervals containing poly-A stretches from learning."));
    addOption(parser, ArgParseOption("ea2", "epaa", "Exclude intervals containing poly-A stretches from analysis."));
    addOption(parser, ArgParseOption("et1", "eptl", "Exclude intervals containing poly-U stretches from learning."));
    addOption(parser, ArgParseOption("et2", "epta", "Exclude intervals containing poly-U stretches from analysis."));
 
    addOption(parser, ArgParseOption("mrtf", "mrtf", "Fit gamma shape k only for positions with min. covariate value.", ArgParseArgument::DOUBLE));
    addOption(parser, ArgParseOption("mtc", "mtc", "Maximum number of truncations at one position. For sites with counts above threshold the whole interval will be discarded! Default: 250.", ArgParseArgument::INTEGER));
    setMinValue(parser, "mtc", "50");
    setMaxValue(parser, "mtc", "254"); 

    addOption(parser, ArgParseOption("pet", "pet", "Prior enrichment threshold: a KDE threshold corresponding to 7 read start counts at one position will be used for initial classification of 'non-enriched' and 'enriched' site. Default: 7", ArgParseArgument::INTEGER));
    setMinValue(parser, "pet", "2");
    setMaxValue(parser, "pet", "50");

    addSection(parser, "General user options");
    addOption(parser, ArgParseOption("nt", "nt", "Number of threads used for learning.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("nta", "nta", "Number of threads used for applying learned parameters. Increases memory usage, if greater than number of chromosomes used for learning, since HMM will be build for multiple chromosomes in parallel.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("tmp", "tmp", "Path to directory to store intermediate files. Default: /tmp ?", ArgParseArgument::STRING));

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));


    //addOption(parser, ArgParseOption("enk", "enk", "Estimate binomial N from KDEs values."));
    addOption(parser, ArgParseOption("nsx", "nsx", "Do not use GSCL simplex2 to estimate Gamma parameters theta and k. Instead update parameters sequentially using the Brents methods. "));
    //addOption(parser, ArgParseOption("ulr", "ulr", "Use log RPKM values as input."));
    //addOption(parser, ArgParseOption("dis", "dis", "Discard intervals containing only one read start."));
    //addOption(parser, ArgParseOption("bs", "bins", "Binsize: currently used to output kde - bin truncCount relationship.", ArgParseArgument::INTEGER));
    //setMinValue(parser, "bins", "1");
    //setMaxValue(parser, "bins", "1000"); 


    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBpureclip\\fP \\fB-i target.bam\\fP \\fB-bai target.bai\\fP \\fB-g ref.fasta\\fP \\fB-o called_crosslinksites.bed\\fP \\fB-nt 10\\fP  \\fB-iv '1;2;3;'\\fP",
                "");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.bamFileName, parser, "in");
    getOptionValue(options.baiFileName, parser, "bai");
    getOptionValue(options.refFileName, parser, "genome");
    getOptionValue(options.outFileName, parser, "out");
    getOptionValue(options.outRegionsFileName, parser, "or");
    getOptionValue(options.parFileName, parser, "par");
    getOptionValue(options.rpkmFileName, parser, "is");
    getOptionValue(options.inputBamFileName, parser, "ibam");
    getOptionValue(options.inputBaiFileName, parser, "ibai");
    if ((options.rpkmFileName != "" && options.inputBamFileName != "") || 
            (options.rpkmFileName != "" && options.inputBaiFileName != "") || 
            (options.inputBamFileName != "" && options.inputBaiFileName == "") ||
            (options.inputBamFileName == "" && options.inputBaiFileName != "") )
    {
        std::cout << "ERROR: If using background signal as covariates, either -is or -ibam and -ibai must be given!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    if (options.rpkmFileName != "" || options.inputBamFileName != "")
        options.useCov_RPKM = true;
    getOptionValue(options.fimoFileName, parser, "fis");
    if (options.fimoFileName != "")
        options.useFimoScore = true;

    getOptionValue(options.intervals_str, parser, "inter");
    if (isSet(parser, "vtb"))
        options.posteriorDecoding = false;
    getOptionValue(options.maxIter_brent, parser, "mibr");
    getOptionValue(options.maxIter_bw, parser, "mibw");
    getOptionValue(options.g1_kMin, parser, "g1kmin");
    getOptionValue(options.g1_kMax, parser, "g1kmax");
    getOptionValue(options.g2_kMin, parser, "g2kmin");
    getOptionValue(options.g2_kMax, parser, "g2kmax");
    //if (isSet(parser, "g1g2k"))
    //    options.g1_k_le_g2_k = true;
    getOptionValue(options.bandwidth, parser, "bdw");

    getOptionValue(options.useKdeThreshold, parser, "mkde");

    getOptionValue(options.nThresholdForP, parser, "ntp");
    getOptionValue(options.minTransProbCS, parser, "mtp");
    getOptionValue(options.maxkNratio, parser, "mkn");
    getOptionValue(options.distMerge, parser, "dm");

    getOptionValue(options.polyAThreshold, parser, "pat");
    if (isSet(parser, "epal"))
        options.excludePolyAFromLearning = true;
    if (isSet(parser, "epaa"))
        options.excludePolyA = true;
    if (isSet(parser, "eptl"))
        options.excludePolyTFromLearning = true;
    if (isSet(parser, "epta"))
        options.excludePolyT = true;
 
    getOptionValue(options.minRPKMtoFit, parser, "mrtf");
    if (isSet(parser, "mrtf"))
        options.mrtf_kdeSglt = false;

    getOptionValue(options.maxTruncCount, parser, "mtc");
 
    getOptionValue(options.nInputMotifs, parser, "nim");

    getOptionValue(options.prior_enrichmentThreshold, parser, "pet");

    getOptionValue(options.numThreads, parser, "nt");
    getOptionValue(options.numThreadsA, parser, "nta");
    getOptionValue(options.tempPath, parser, "tmp");

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;


    //if (isSet(parser, "enk"))
    //    options.estimateNfromKdes = true;
    if (isSet(parser, "nsx"))
        options.gslSimplex2 = false;   
    //if (isSet(parser, "ulr"))
    //    options.useLogRPKM = true;
    //if (isSet(parser, "dis"))
    //    options.discardSingletonIntervals = true;
    //getOptionValue(options.binSize, parser, "bins");
    getOptionValue(options.applyChr_str, parser, "chr");

    return ArgumentParser::PARSE_OK;
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

    std::cout << "Protein-RNA crosslink site detection \n"
              << "===============\n\n";

    // Print the command line arguments back to the user.
    /*if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << "\n\n";
    }*/
    /////////////////////////////////////////

    if (options.useCov_RPKM)
    {
        GAMMA2_REG gamma1; 
        GAMMA2_REG gamma2;

        if (options.useFimoScore)
        {
            ZTBIN_REG bin1;
            ZTBIN_REG bin2;
            bin1.b0 = log(options.p1/(1.0 - options.p1));            
            bin2.b0 = log(options.p2/(1.0 - options.p2)); 
            resize(bin1.regCoeffs, options.nInputMotifs, 0.0, Exact());
            resize(bin2.regCoeffs, options.nInputMotifs, 0.0, Exact());
            doIt(gamma1, gamma2, bin1, bin2, options);  
        }
        else
        {
            ZTBIN bin1;
            ZTBIN bin2;
            bin1.p = options.p1;            
            bin2.p = options.p2; 
            doIt(gamma1, gamma2, bin1, bin2, options);  
        }
    }
    else
    {
        GAMMA2 gamma1;           
        GAMMA2 gamma2;

        if (options.useFimoScore)
        {
            ZTBIN_REG bin1;
            ZTBIN_REG bin2;

            bin1.b0 = log(options.p1/(1.0 - options.p1));            
            bin2.b0 = log(options.p2/(1.0 - options.p2)); 
            resize(bin1.regCoeffs, options.nInputMotifs, 0.0, Exact());
            resize(bin2.regCoeffs, options.nInputMotifs, 0.0, Exact());
            doIt(gamma1, gamma2, bin1, bin2, options);  
        }
        else
        {
            ZTBIN bin1;
            ZTBIN bin2;
            bin1.p = options.p1;            
            bin2.p = options.p2; 
            doIt(gamma1, gamma2, bin1, bin2, options);  
        }

    }
    return 0;
}






