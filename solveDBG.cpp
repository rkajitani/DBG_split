/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus-allee.

Platanus-allee is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus-allee is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus-allee; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "solveDBG.h"
#include "seqlib.h"
#include "kmer.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>
#include <cfloat>
#include <iomanip>

using std::vector;
using std::string;
using std::unordered_map;
using std::cerr;
using std::endl;


//////////////////////////////////////////////////////////////////////////////////////
// const parameter define
//////////////////////////////////////////////////////////////////////////////////////
const double SolveDBG::MIN_LONG_READ_LENGTH_CUTOFF_FACTOR = 1;
const double SolveDBG::MAX_LONG_READ_LENGTH_CUTOFF_FACTOR = 4;
const double SolveDBG::LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR = 0.1;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_COVERAGE = 0.8;

//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
SolveDBG::SolveDBG()
: Scaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-k"] = "1";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-i"] = "0.8";
    optionSingleArgs["-unlink_list"] = "";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-anchor_bubble"] = vector<string>();
    optionMultiArgs["-anchor_homo"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-ont"] = vector<string>();
    optionMultiArgs["-gc"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    optionMultiArgs["-S"] = vector<string>(2);
    optionMultiArgs["-S"][0] = "500";
    optionMultiArgs["-S"][1] = "1000";

    optionBool["-no_scaffold"] = false;
    optionBool["-unphase"] = false;
    optionBool["-reduce_redundancy"] = false;
	optionBool["-divide_only"] = false;
	optionBool["-simplify_junction"] = false;

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";

    optionSingleArgs["-mapper"] = "";
    optionBool["-minimap2_sensitive"] = false;
    optionBool["-minimap2"] = true;
    optionBool["-minialign"] = false;
    optionBool["-kmer_align"] = false;
    optionBool["-emem"] = false;

    optionBool["-fastg"] = false;
    optionBool["-aggressive"] = false;
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::usage(void) const
{
    std::cerr << "\nUsage: platanus_allee solveDBG [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n"
              << "    -anchor_bubble FILE1 [FILE2 ...]   : anchor bubble_seq_file (fasta format)\n"
              << "    -anchor_homo FILE1 [FILE2 ...]     : anchor homo_seq_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)\n"
              << "    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)\n"
              << "    -gc FILE1 [FILE2 ...]              : Guiding contig file; i.e. other assemblies, synthetic long-reads or corrected reads (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : tagged_pair_files (10x Genomics) (reads in 1 file, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : tagged_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)\n"
              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
              << "    -a{INT} INT                        : lib_id average_insert_size\n"
              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -e FLOAT                           : coverage depth of homozygous region (default auto)\n"
              << "    -L INT                             : maximum fragment length of tag (10x Genomics) (default " << optionSingleArgs.at("-L") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1] << " " << optionMultiArgs.at("-s")[2]  << ")\n"
              << "    -S INT1 [INT2 ...]                 : minimum alignment length for long reads (default " << optionMultiArgs.at("-S")[0] << ")\n"
              << "    -i FLOAT                           : minimum identity for long-read alignment (identity, default " << optionSingleArgs.at("-i") << ")\n"
              << "    -k INT                             : minimum number of links to phase variants (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -unphase                           : not phase heterozygous regions and construct consensus scaffolds (default false)\n"
              << "    -simplify_junction                 : only simplify junction-contigs based on DBG-overlaps (default false)\n"
              << "    -mapper FILE                       : path of mapper executable file (default, minimap2; only effective with -p option)\n"
              << "    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)\n"
              << "    -reduce_redundancy                 : reduce redundant sequences that exactly matche others (default, off)\n"
              << "    -divide_only                       : only divide input sequences (default, off)\n"
              << "    -unlink_list                       : reduce redundant sequences that exactly matche others (default, off)\n"
//              << "    -minimap2                          : use minimap2 insterd of minimap (default) for alignment of long reads\n"
//              << "    -minialign                         : use minialign insterd of minimap (default) for alignment of long reads\n"
//              << "    -kmer_align                        : use built-in alignment based on k-mer match insterd of minimap (default) for long reads\n"
//              << "    -fastg                             : output only fastg files of graphs for Bandage (default off)\n"
              << "    -emem                              : use E-mem insterd of minimap2 for alignment of long reads\n"
              << "    -aggressive                        : aggressively extend scaffolds (default off)\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x and -X options.\n"
              << "\n\n"


              << "Outputs:\n"
			  << "    PREFIX_*.fa\n"
              << "\n"
              << "Outputs (-fastg):\n"
			  << "    PREFIX_lib#_graph.fastg\n"
             << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initialize parameters
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::initializeParameters(void)
{
	contigMaxK = platanus::Contig::getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
	seedLength = std::min((unsigned)contigMaxK - 1, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);

    multiSeedLengthForShortRead.clear();
	if (!(optionPairFile.empty())) {
		for (auto itr = optionMultiArgs["-s"].begin(); itr != optionMultiArgs["-s"].end(); ++itr) {
			unsigned length = std::stoi(*itr);
			multiSeedLengthForShortRead.push_back(length);
			if (seedLength > length)
				seedLength = length;
		}
	}
	if ((optionBool["-emem"] || optionBool["-kmer_align"]) && !(optionMultiArgs["-p"].empty())) {
		for (auto itr = optionMultiArgs["-S"].begin(); itr != optionMultiArgs["-S"].end(); ++itr) {
			unsigned length = std::stoi(*itr);
			multiSeedLengthForLongRead.push_back(length);
			if (seedLength > length && optionBool["-kmer_align"])
				seedLength = length;
		}
	}

    minAlignmentLongRead.clear();
	for (auto itr = optionMultiArgs["-S"].begin(); itr != optionMultiArgs["-S"].end(); ++itr) {
		unsigned length = std::stoi(*itr);
		minAlignmentLongRead.push_back(length);
	}
    sort(minAlignmentLongRead.begin(), minAlignmentLongRead.end());

    keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);
    bubbleThreshold = atof(optionSingleArgs["-u"].c_str());
    minIdentity = atof(optionSingleArgs["-i"].c_str());
    minLink = atoi(optionSingleArgs["-l"].c_str());
    minLinkToPhase = atoi(optionSingleArgs["-k"].c_str());
    minOverlapForScaffolding = atoi(optionSingleArgs["-v"].c_str());
    numThread = atoi(optionSingleArgs["-t"].c_str());
    pairedDBG.setSeedLength(seedLength);
    pairedDBG.setMinTolerenceFactor(MIN_TOL_FACTOR);
    pairedDBG.setMaxFragmentLengthOfTag(atoi(optionSingleArgs["-L"].c_str()));

    unsigned numLibrary = 0;
	if (optionMultiArgs["-anchor_bubble"].empty() && optionMultiArgs["-anchor_homo"].empty()) {
		sort(optionPairFile.begin(), optionPairFile.end());
		numFilePerLibraryID.resize(optionPairFile.size());
		libraryIDList.resize(optionPairFile.size());
		for (unsigned i = 0; i < optionPairFile.size(); ++i) {
			++(numFilePerLibraryID[numLibrary]);
			libraryIDList[numLibrary] = optionPairFile[i].libraryID;
			if (i + 1 >= optionPairFile.size() || optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
				++numLibrary;
			}
		}
	}
	else {
		numLibrary = 1;
	}
    libraryMT.resize(numLibrary);

    omp_set_num_threads(numThread);
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}


//////////////////////////////////////////////////////////////////////////////////////
// exec scaffold
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::exec(void)
{
    initializeParameters();

	if (!(optionMultiArgs["-anchor_bubble"].empty() && optionMultiArgs["-anchor_homo"].empty())) {
		PairedDBG dividedDBG;
		mapLibraryAndInitGraphAnchorMode(dividedDBG, numThread);
		dividedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		dividedDBG.makeGraph(numThread);
		dividedDBG.adjustOppositeBubbleNodeIDDirection();
		dividedDBG.setTolerence(this->contigMaxK);
		dividedDBG.loadResultSeqSimple(this->contigMaxK, this->contigReadLength, "seq");
		dividedDBG.markRedundantResultSeq(numThread);
		dividedDBG.outputResultSeqWithBubble(optionSingleArgs["-o"], "_anchoredPrimaryBubble.fa", "_anchoredSecondaryBubble.fa", "_anchoredPrimaryFork.fa", "_anchoredSecondaryFork.fa", "_anchoredNestedBubble.fa", "_unanchoredNonBubbleOther.fa", "_anchoredBubbleRelation.tsv", this->contigMaxK, this->contigReadLength);

		std::ostringstream oss;
		oss << "rm " << optionSingleArgs["-o"] << "_anchoredPrimaryFork.fa "
			<< optionSingleArgs["-o"] << "_anchoredSecondaryFork.fa "
			<< "; "
			<< "cat " << optionSingleArgs["-o"] << "_anchoredPrimaryBubble.fa "
			<< optionSingleArgs["-o"] << "_anchoredSecondaryBubble.fa "
			<< optionSingleArgs["-o"] << "_anchoredNestedBubble.fa "
			<< optionSingleArgs["-o"] << "_unanchoredNonBubbleOther.fa "
			<< ">" << optionSingleArgs["-o"] <<  "_AllAnchorResults.fa"
		;
		system(oss.str().c_str());

		cerr << "solve_DBG completed!" << endl;
		return;
	}

    mapLibraryAndInitGraph(numThread);

	if (optionBool["-fastg"]) {
		generateGraphFastg();
		cerr << "scaffold completed!" << endl;
		return;
	}

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.makeGraph(numThread);


	if (optionBool["-simplify_junction"]) {
		pairedDBG.joinUnambiguousNodePairIterative(numThread);

		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "seq");
		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_simplified.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_simplifiedComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	if (optionBool["-unphase"]) {
		if (optionSingleArgs["-e"] == "")
			pairedDBG.calculateHeteroAndAverageCoverageUnphase();
		else
			pairedDBG.setHeteroCoverage(atof(optionSingleArgs["-e"].c_str()));

		pairedDBG.clearEdges();

		extendConsensus(4, !(optionBool["-aggressive"]));

		if (!(libraryMT.empty()))
			pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
		else
			pairedDBG.setTolerence(this->contigMaxK);
		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");

		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_consensusScaffold.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_consensusScaffoldComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	if (optionBool["-divide_only"]) {
		if (optionSingleArgs["-e"] == "")
			pairedDBG.calculateHeteroAndAverageCoverageUnphase();
		else
			pairedDBG.setHeteroCoverage(atof(optionSingleArgs["-e"].c_str()));

		pairedDBG.clearEdges();

		extendConsensusToEstimateInsertSize();
		pairedDBG.resetGraph();

		pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, true, numThread);

		pairedDBG.loadDividedContigResultSeq(this->contigMaxK, this->contigReadLength);
		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_divided.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_dividedComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	pairedDBG.extractDBGBubbleInformation();
	pairedDBG.clearEdges();

	extendConsensusToEstimateInsertSize();
	pairedDBG.resetGraph();


	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setCutoffLength(0);
	pairedDBG.makeGraph(numThread);

	pairedDBG.extractDBGBubbleInformation();
//	pairedDBG.setOppositeForkContigIDOverlapped(numThread);
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.setOppositeBubbleContigIDByEndMatch();
//	pairedDBG.setOppositeBubbleContigIDByOneEndMatch();

//	pairedDBG.setForkJunctionContigIDOverlapped();
	pairedDBG.setBubbleJunctionContigIDOverlapped();

	pairedDBG.clearEdges();


	for (long outerIteration = 0; outerIteration < 4; ++outerIteration) {
		if (longReadLibraryMT.size() > 0)
			pairedDBG.loadLongReadLink(minAlignmentLongRead.front(), LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, this->contigMaxK, numThread);

		for (long iteration = 0; iteration < 2; ++iteration) {
			if (optionBool["-no_scaffold"])
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::BUBBLE_AWARE_MODE);
			else if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			pairedDBG.setCutoffLength(0);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setTargetLibraryIndex(i);
				unsigned tolerenceFactor = MAX_TOL_FACTOR;
				pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
				cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
			}

			if (longReadLibraryMT.size() > 0) {
				if (optionBool["-no_scaffold"])
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::BUBBLE_AWARE_MODE);
				else if (iteration == 0)
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				else
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize() << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			if (optionBool["-no_scaffold"])
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::BUBBLE_AWARE_MODE);
			else if (iteration == 0 || optionBool["-no_scaffold"])
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			if (optionBool["-no_scaffold"]) {
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
				pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			}
		}

		if (optionBool["-no_scaffold"]) {
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			continue;
		}


		pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
		pairedDBG.setCutoffLength(0);
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);
		pairedDBG.clearEdges();


		for (long iteration = 0; iteration < 2; ++iteration) {
			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();
					if (iteration > 0)
						pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}

			if (longReadLibraryMT.size() > 0) {
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageLength() << endl;
				pairedDBG.setTolerence(2 * this->contigMaxK);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR * longReadLibraryMT[0].getAverageInsSize());
					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}

			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					if (iteration == 0)
						pairedDBG.setMinLink(minLinkToPhase);
					else {
						pairedDBG.setMinLink(minLink);
						pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					}
					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}
		}


		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleContigInNonHeteroNode();
		pairedDBG.divideBubbleJunctionNode(false);


		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				long linkThreshold;
				if (iteration % 2 == 0)
					linkThreshold = minLink;
				else
					linkThreshold = std::max(minLink, pairedDBG.estimateLink());

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();


					pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					pairedDBG.setMinLink(minLink);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();
				}
			}
			if (outerIteration < 2)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(false, numThread);


		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMinOverlap(minOverlapForScaffolding);

			for (long iteration = 0; iteration < 2; ++iteration) {
				long alignmentIndex = minAlignmentLongRead.size() - std::min(iteration + 1, (long)minAlignmentLongRead.size());
				pairedDBG.loadLongReadLink(minAlignmentLongRead[alignmentIndex], LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, this->contigMaxK, numThread);

				pairedDBG.divideNodeUsingBubbleContigPair(numThread);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					if (outerIteration == 0)
						pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
					else
						pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

					pairedDBG.setTolerence(2 * this->contigMaxK);
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());

					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();

					pairedDBG.setMinLink(minLink);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffoldCombine();

					pairedDBG.setMinLink(minLink);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

					pairedDBG.setMinLink(minLink);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffoldCombine();

					pairedDBG.divideNodeUsingBubbleContigPair(numThread);
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
					pairedDBG.setMinLink(minLinkToPhase);
					pairedDBG.setCutoffLength(0);
					pairedDBG.joinUnambiguousNodePairIterative(numThread);
					pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}

			pairedDBG.setMinOverlap(this->contigMaxK - 1);
		}


		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			if (longReadLibraryMT.size() > 0) {
				long alignmentIndex = minAlignmentLongRead.size() - std::min(iteration + 1, (long)minAlignmentLongRead.size());
				pairedDBG.loadLongReadLink(minAlignmentLongRead[alignmentIndex], LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, this->contigMaxK, numThread);
			}

			long linkThreashold = (iteration + 1) * minLink;
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);

					pairedDBG.trimSparseEnd();
					pairedDBG.trimRepeatEnd();

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();

					pairedDBG.setMinLink(minLink);
					pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteConflictingBubbleEdge(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();

//					if (outerIteration < 2)
//						pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
				}
			}
			if (outerIteration < 2)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleJunctionNode(true);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);

		if (outerIteration < 2) {
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
		}
	}


	if (!optionBool["-aggressive"])
		pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, true, true, false, numThread);

	pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);


	pairedDBG.divideBubbleContigInNonHeteroNode();

	pairedDBG.setMode(0);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	pairedDBG.copyAllNodes(phasedGraph);

	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();


	extendConsensus(1, false);

	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");
	outputGraph("_preliminaryConsensusScaffold.fa");
	pairedDBG.clearResultSeq();

	pairedDBG.setMode(0);
	pairedDBG.remakeGraphRecoveringSecondaryBubble(phasedGraph);
	pairedDBG.makeGraph(numThread);
	pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);

    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");

	if (optionBool["-reduce_redundancy"])
		pairedDBG.markRedundantResultSeq(numThread);

	pairedDBG.outputResultSeqWithBubble(optionSingleArgs["-o"], "_primaryBubble.fa", "_secondaryBubble.fa", "_primaryFork.fa", "_secondaryFork.fa", "_nestedBubble.fa", "_nonBubbleOther.fa", "_bubbleRelation.tsv", this->contigMaxK, this->contigReadLength);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");
	pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_phasedScaffoldComponent.bed");
	
	cerr << "solve_DBG completed!" << endl;
	return;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::mapLibraryAndInitGraph(const int numThread)
{
    platanus::Contig contig;
    std::unique_ptr<HeteroMapper> mapper(new HeteroMapper(seedLength, keyLength));

    readLibrary(mapper, contig, numThread);
    cerr << "CONTIG_AVERAGE_COVERAGE = " << averageCoverage << endl;

	if (optionBool["-fastg"]) {
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_nameIndex.csv";
		contig.printNameIndexCSV(outStream.str());
	}

	mapper->setMultiSeedLength(multiSeedLengthForShortRead);
    unsigned nowFileNumber = 0;
    for (unsigned i = 0; i < libraryMT.size(); ++i) {
		libraryMT[i][0].setAverageInsSize(0);
        int nowLibraryID = optionPairFile[nowFileNumber].libraryID;
        cerr << "[LIBRARY " << libraryIDList[i] << "]" << endl;
        // set average length and minimum insert size
        // estimate insert size
        libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));
        long minInsertion = optionMinIns.find(nowLibraryID) == optionMinIns.end() ? 0 : optionMinIns[nowLibraryID];

        mapper->contigMap.mapPairAndSaveReadLink(libraryMT[i], minInsertion, this->contigMaxK + 1, numThread);

		for (int j = 0; j < numThread; ++j) {
			fclose(libraryMT[i][j].pairFP);
			libraryMT[i][j].pairFP = NULL;
		}

        libraryMT[i][0].setInsCutoffRate(optionInsCutoffRate.find(nowLibraryID) == optionInsCutoffRate.end() ? DEFAULT_INS_CUTOFF_RATE : optionInsCutoffRate[nowLibraryID]);
        if (optionAveIns.find(nowLibraryID) != optionAveIns.end() || optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
            if (optionAveIns.find(nowLibraryID) != optionAveIns.end()) {
                libraryMT[i][0].setAverageInsSize(optionAveIns[nowLibraryID]);
                std::cerr << "Average insert size specified: AVE = " << libraryMT[i][0].getAverageInsSize() << std::endl;
            }
            if (optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
                libraryMT[i][0].setSDInsSize(optionSDIns[nowLibraryID]);
            } else {
                libraryMT[i][0].setSDInsSize(static_cast<long>(static_cast<double>(libraryMT[i][0].getAverageInsSize()) / 10.0 + 0.5));
            }
        }
        nowFileNumber += numFilePerLibraryID[i];
    }

	if (longReadLibraryMT.size() > 0) {
		cerr << "[LONG_READ LIBRARY]" << endl;

		string alignerOutFilename(optionSingleArgs["-o"]);
		alignerOutFilename += "_longReadAlignment.tsv";

		string aligner = optionSingleArgs["-mapper"];

		if (optionBool["-kmer_align"]) {
			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			mapper->setMultiSeedLength(multiSeedLengthForLongRead);
			mapper->contigMap.mapLongReadAndSaveReadLink(longReadLibraryMT, this->contigMaxK, numThread);
		}
		else if (optionBool["-minialign"]) {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "minialign";

			string alignerTargetFilename(optionSingleArgs["-o"]);
			alignerTargetFilename += "_longReadTarget.fa";
			std::ostringstream targetOss;
			targetOss << "cat ";
			for (auto itr = optionMultiArgs["-c"].begin(); itr != optionMultiArgs["-c"].end(); ++itr)
				targetOss << " " << *itr;
			for (auto itr = optionMultiArgs["-b"].begin(); itr != optionMultiArgs["-b"].end(); ++itr)
				targetOss << " " << *itr;
			targetOss << " >" << alignerTargetFilename;
			system(targetOss.str().c_str());

			string alignerQueryFilename(optionSingleArgs["-o"]);
			alignerQueryFilename += "_longReadQuery.txt";
			std::ostringstream queryOss;
			queryOss << "cat ";
			for (auto itr = optionMultiArgs["-p"].begin(); itr != optionMultiArgs["-p"].end(); ++itr)
				queryOss << " " << *itr;
			queryOss << " >" << alignerQueryFilename;
			system(queryOss.str().c_str());

			execMinialign(alignerTargetFilename, alignerQueryFilename, alignerOutFilename , numThread, aligner);

			std::ostringstream rmOss;
			rmOss << "rm " << alignerTargetFilename << " " << alignerQueryFilename;
			system(rmOss.str().c_str());

			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			readLongReadPAFfileAndSaveLink(alignerOutFilename, contig, longReadLibraryMT, minAlignmentLongRead.front(), LONG_READ_MIN_ALIGNMENT_COVERAGE, 0, contigMaxK, numThread);
		}
		else if (optionBool["-emem"]) {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "e-mem-paf";
			execEmem(optionMultiArgs["-c"], optionMultiArgs["-b"], optionMultiArgs["-p"], alignerOutFilename, multiSeedLengthForLongRead.front(), numThread, aligner);
			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			readLongReadPAFfileAndSaveLink(alignerOutFilename, contig, longReadLibraryMT, multiSeedLengthForLongRead.front(), 0.0, 0.0, contigMaxK, numThread);
		}
		else if (optionBool["-minimap2"]) {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "minimap2";

			vector<string> targetFilename = optionMultiArgs["-c"];
			execMinimap2(targetFilename, optionMultiArgs["-b"], pacBioLongReadFilename, alignerOutFilename, numThread, aligner, "-x map-hifi --secondary=no");
			execMinimap2(targetFilename, optionMultiArgs["-b"], nanoporeLongReadFilename, alignerOutFilename, numThread, aligner, "-x map-ont --secondary=no");
			execMinimap2(targetFilename, optionMultiArgs["-b"], guideContigLongReadFilename, alignerOutFilename, numThread, aligner, "-x asm10 --secondary=no");

			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			readLongReadPAFfileAndSaveLink(alignerOutFilename, contig, longReadLibraryMT, minAlignmentLongRead.front(), LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, contigMaxK, numThread);
		}
		else {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "minimap";
			execMinimap(optionMultiArgs["-c"], optionMultiArgs["-b"], optionMultiArgs["-p"], alignerOutFilename, contigMaxK/2 , numThread, aligner);
			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			readLongReadPAFfileAndSaveLink(alignerOutFilename, contig, longReadLibraryMT, minAlignmentLongRead.front(), LONG_READ_MIN_ALIGNMENT_COVERAGE, 0.0, contigMaxK, numThread);
		}

		for (int j = 0; j < numThread; ++j) {
			fclose(longReadLibraryMT[j].pairFP);
			longReadLibraryMT[j].pairFP = NULL;
		}

		vector<long> insSizeDistribution;
		longReadLibraryMT[0].readInsertSizeFile(insSizeDistribution);
		longReadLibraryMT[0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_longReadLibrary" << "_readDistribution.tsv";
		longReadLibraryMT[0].printInsertSizeFreq(outStream.str());
		cerr << "[LONG_READ_LIBRARY " << 1 << "]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize()
			 << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
	}	

	if (tagLibraryMT.size() > 0) {
		mapper->setMultiSeedLength(multiSeedLengthForShortRead);

		cerr << "[TAG LIBRARY]" << endl;
		mapper->contigMap.mapTagPairMT(tagLibraryMT, numThread);

		for (int j = 0; j < numThread; ++j) {
			fclose(tagLibraryMT[j].pairFP);
			tagLibraryMT[j].pairFP = NULL;
		}
	}	

	if (libraryMT.size() > 0) {
		pairedDBG.setAllLibraryMT(&libraryMT);
		pairedDBG.setTargetLibraryIndex(0);
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);
		pairedDBG.setContigNameIndex(contig.nameIndex);

		if (optionSingleArgs["-unlink_list"].size() > 0)
			pairedDBG.setContigUnlinkFlags(optionSingleArgs["-unlink_list"]);

		vector<long> insSizeDistribution;
		libraryMT[0][0].readInsertSizeFile(insSizeDistribution);
		pairedDBG.insertSizeDistribution(libraryMT[0], insSizeDistribution, numThread);

		vector<long> seqLengths;
		pairedDBG.scaffoldLengthList(seqLengths);

		if (libraryMT[0][0].getAverageInsSize() == 0)
			libraryMT[0][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << 1 << "_insFreq.tsv";
		libraryMT[0][0].printInsertSizeFreq(outStream.str());
		cerr << "[LIBRARY " << 1 << "]\nAVE_INS = " << libraryMT[0][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[0][0].getSDInsSize() << endl;
	}
	else {
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);
	}

    pairedDBG.setContigMaxK(this->contigMaxK);
    pairedDBG.setMinOverlap(this->contigMaxK - 1);
    pairedDBG.saveOverlap(mapper->contigMap, this->contigMaxK - 1, this->contigMaxK, numThread);
    pairedDBG.classifyNode();

	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setLongReadLibraryMT(&longReadLibraryMT);
		pairedDBG.loadLongReadLink(minAlignmentLongRead.front(), LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, this->contigMaxK, numThread);
	}

	if (tagLibraryMT.size() > 0) {
		pairedDBG.setTagLibraryMT(&tagLibraryMT);
		pairedDBG.countMappedTagForEachContig(numThread);
	}

    cerr << "destructing mapper objects..." << std::endl;
}


void SolveDBG::mapLibraryAndInitGraphAnchorMode(PairedDBG &dividedDBG, const int numThread)
{
    platanus::Contig contig;
    platanus::Contig anchorBubble;
    platanus::Contig anchorHomo;
	vector<platanus::Region> anchorBubbleMap;
	vector<platanus::Region> anchorHomoMap;
    std::unique_ptr<HeteroMapper> mapper(new HeteroMapper(seedLength, keyLength));

    readLibraryAnchorMode(mapper, contig, anchorBubble, anchorHomo, numThread);
	mapper->contigMap.setSeedLength(this->contigMaxK);
    cerr << "ANCHOR_CONTIG_AVERAGE_COVERAGE = " << averageCoverage << endl;

	if (!optionMultiArgs["-anchor_bubble"].empty()) {
		cerr << "[LIBRARY ANCHOR_BUBBLE]" << endl;
		mapper->contigMap.mapAnchorBubbleMT(anchorBubble, anchorBubbleMap, numThread);
	}

	if (!optionMultiArgs["-anchor_homo"].empty()) {
		cerr << "[LIBRARY ANCHOR_HOMO]" << endl;
		mapper->contigMap.mapAnchorHomoMT(anchorHomo, anchorHomoMap, numThread);
	}

	pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
	pairedDBG.setContigName(contig.name);
    pairedDBG.setContigMaxK(this->contigMaxK);
    pairedDBG.setMinOverlap(this->contigMaxK - 1);

	dividedDBG.reconstructAnchorDividedGraph(pairedDBG, contig, anchorBubble, anchorHomo, anchorBubbleMap, anchorHomoMap);

    cerr << "destructing mapper objects..." << std::endl;
}


void SolveDBG::readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread)
{
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1)
    for (int i = -3; i < static_cast<int>(libraryMT.size()); ++i) {
        try {
            if (i == -3) {
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-c"][i]);

				long j = contig.numSeq;
                for (unsigned i = 0; i < optionMultiArgs["-b"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-b"][i]);

				pairedDBG.setNumInputBubbleContig(contig.numSeq - j);
				contig.setNameIndex();

                this->contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
                this->contigReadLength = contig.getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
                if (this->contigReadLength == 0)
                    this->contigReadLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;

				if (optionSingleArgs["-e"] == "")
					averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
				else
					averageCoverage = atof(optionSingleArgs["-e"].c_str());

                mapper->setContigMap(contig);
                mapper->makeKmerTableContigMap();
            }
			else if (i == -2) {
                if (optionMultiArgs["-p"].empty() && optionMultiArgs["-ont"].empty() && optionMultiArgs["-gc"].empty())
					continue;

				vector<string> filenames(optionMultiArgs["-p"]);
				pacBioLongReadFilename = optionMultiArgs["-p"];

				std::copy(optionMultiArgs["-ont"].begin(), optionMultiArgs["-ont"].end(), std::back_inserter(filenames));
				nanoporeLongReadFilename = optionMultiArgs["-ont"];

				std::copy(optionMultiArgs["-gc"].begin(), optionMultiArgs["-gc"].end(), std::back_inserter(filenames));
				guideContigLongReadFilename = optionMultiArgs["-gc"];

				longReadLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					longReadLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(filenames.size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(filenames[j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleMT(longReadLibraryMT, filenames[j], numThread, false, isFastq, true);
				}
            }
			else if (i == -1) {
                if (optionMultiArgs["-x"].size() == 0 && optionMultiArgs["-X"].size() == 0)
					continue;

				vector<string> filenames(optionMultiArgs["-x"]);
				filenames.insert(filenames.end(), optionMultiArgs["-X"].begin(), optionMultiArgs["-X"].end());

				std::unordered_map<string, int> tagStringConverter;
				setTagStringConverter(filenames, tagStringConverter);

				tagLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					tagLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-x"].size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-x"][j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleTaggedMT(tagLibraryMT, optionMultiArgs["-x"][j], numThread, false, isFastq, true, tagStringConverter);
				}

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-X"].size()); j += 2) {
					platanus::FILETYPE fileFormat1 = checkFileFormat(optionMultiArgs["-X"][j]);
					platanus::FILETYPE fileFormat2 = checkFileFormat(optionMultiArgs["-X"][j + 1]);
					if (fileFormat1 != fileFormat2) {
						throw platanus::FormatError("Different file type in paired-file (-X).");
					}
					bool isFastq;
					switch (fileFormat1) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaPairTaggedMT(tagLibraryMT, optionMultiArgs["-X"][j] , optionMultiArgs["-X"][j + 1], numThread, false, isFastq, tagStringConverter);
				}
            }
			else if (!optionPairFile.empty()) {
                unsigned nowFileNumber = 0;
                libraryMT[i].resize(numThread);
                for (int j = 0; j < i; ++j)
                    nowFileNumber += numFilePerLibraryID[j];

                for (int j = 0; j < numThread; ++j) {
                    libraryMT[i][j].pairFP = platanus::makeTemporaryFile();
                }
                for (int j = 0; j < numFilePerLibraryID[i]; ++j) {
                    bool isFastq, isMate;
                    isMate = optionPairFile[nowFileNumber+j].libraryType == "-op"
                           || optionPairFile[nowFileNumber+j].libraryType == "-OP" ? true : false;
                    if (optionPairFile[nowFileNumber+j].libraryType == "-ip" || optionPairFile[nowFileNumber+j].libraryType == "-op") {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        switch (fileFormat) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaSingleMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, numThread, isMate, isFastq);
                    }
					else {
                        platanus::FILETYPE fileFormat1 = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        platanus::FILETYPE fileFormat2 = checkFileFormat(optionPairFile[nowFileNumber+j].fileSecond);
                        if (fileFormat1 != fileFormat2) {
                            throw platanus::FormatError("Different file type in paired-file.");
                        }
                        switch (fileFormat1) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaPairMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, optionPairFile[nowFileNumber+j].fileSecond, numThread, isMate, isFastq);
                    }
                }
            }
        } catch (platanus::ErrorBase &e) {
            e.showErrorMessage();
            exit(e.getID());
        }
    }
}


void SolveDBG::readLibraryAnchorMode(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, platanus::Contig &anchorBubble, platanus::Contig &anchorHomo, const int numThread)
{
    omp_set_num_threads(numThread);

	# pragma omp parallel sections
	{
		# pragma omp section
		{
			for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i)
				contig.readFastaCoverage(optionMultiArgs["-c"][i]);

			long j = contig.numSeq;
			for (unsigned i = 0; i < optionMultiArgs["-b"].size(); ++i)
				contig.readFastaCoverage(optionMultiArgs["-b"][i]);

			pairedDBG.setNumInputBubbleContig(contig.numSeq - j);
			contig.setNameIndex();

			mapper->setContigMap(contig);
			mapper->makeKmerTableContigMap();
		}

		# pragma omp section
		{
			if (!optionMultiArgs["-anchor_bubble"].empty()) {
				platanus::Contig rawBubble;
				for (unsigned i = 0; i < optionMultiArgs["-anchor_bubble"].size(); ++i)
					rawBubble.readFastaCoverage(optionMultiArgs["-anchor_bubble"][i]);

				rawBubble.setNameIndex();
				pairedDBG.calculateHeteroCoverageContig(rawBubble);
				pairedDBG.filterAnchorBubble(rawBubble, anchorBubble);
				rawBubble.clear();

				pairedDBG.calculateHeteroCoverageContig(anchorBubble);
				this->contigMaxK = anchorBubble.getMaxKFromFastaHeader(optionMultiArgs["-anchor_bubble"][0]);
				this->contigReadLength = anchorBubble.getReadLengthFromFastaHeader(optionMultiArgs["-anchor_bubble"][0]);
			}

			if (!optionMultiArgs["-anchor_homo"].empty()) {
				platanus::Contig rawHomo;
				for (unsigned i = 0; i < optionMultiArgs["-anchor_homo"].size(); ++i)
					rawHomo.readFastaCoverage(optionMultiArgs["-anchor_homo"][i]);

				rawHomo.setNameIndex();
				if (optionMultiArgs["-anchor_bubble"].empty()) {
					pairedDBG.calculateHeteroCoverageContig(rawHomo, true);
				}
				pairedDBG.filterAnchorHomo(rawHomo, anchorHomo);
				rawHomo.clear();

				if (optionMultiArgs["-anchor_bubble"].empty()) {
					pairedDBG.calculateHeteroCoverageContig(anchorHomo, true);
					this->contigMaxK = anchorHomo.getMaxKFromFastaHeader(optionMultiArgs["-anchor_homo"][0]);
					this->contigReadLength = anchorHomo.getReadLengthFromFastaHeader(optionMultiArgs["-anchor_homo"][0]);
				}
			}
		}
	}

	if (this->contigReadLength == 0)
		this->contigReadLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;
	if (optionSingleArgs["-e"] == "")
		averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
	else
		averageCoverage = atof(optionSingleArgs["-e"].c_str());
}




//////////////////////////////////////////////////////////////////////////////////////
// destroy graph
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::outputAndAfterTreatment(void)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    string componentFilename(outFilename);
    outFilename += "_solvedContig.fa";
    componentFilename += "_solvedContigComponent.tsv";
    pairedDBG.cutAndPrintSeq(this->contigMaxK, this->contigReadLength, outFilename, componentFilename);
}

void SolveDBG::outputGraph(const char *suffix)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    outFilename += suffix;
    pairedDBG.printResultSeq(outFilename);
}

void SolveDBG::updateAndWriteInsertSize(const long libraryIndex)
{
	pairedDBG.updateInsertLengthFP(libraryMT[libraryIndex], numThread);
	vector<long> insSizeDistribution;
	libraryMT[libraryIndex][0].readInsertSizeFile(insSizeDistribution);
	pairedDBG.insertSizeDistribution(libraryMT[libraryIndex], insSizeDistribution, numThread);
	if (libraryIndex > 0)
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, libraryMT[libraryIndex - 1][0].getAverageInsSize(), platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);
	else
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

	std::ostringstream outStream;
	outStream << optionSingleArgs["-o"] << "_lib" << (libraryIndex + 1) << "_insFreq.tsv";
	printInsertSizeFreq(outStream.str(), insSizeDistribution);
}

void SolveDBG::extendConsensus(const long numOuterIteration, const bool divisionFlag)
{
	pairedDBG.clearContigPreviousParentNodeID();

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);

	pairedDBG.setMinOverlap(minOverlapForScaffolding);

	for (long outerIteration = 0; outerIteration < numOuterIteration; ++outerIteration) {
		long numIteration = 2;
		for (long iteration = 0; iteration < numIteration; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (long libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
				pairedDBG.setTargetLibraryIndex(libraryIndex);

				if (libraryMT[libraryIndex][0].getAverageInsSize() <= 0)
					updateAndWriteInsertSize(libraryIndex);

				if (iteration == 0)
					pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					if (iteration > 0)
						pairedDBG.deleteErroneousEdgeIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();

					if (numIteration >= 4 && outerIteration < 2)
						pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
				}
			}

			if (divisionFlag)
				pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
		}

		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();


		if (longReadLibraryMT.size() > 0) {
			for (long iteration = 0; iteration < numIteration; ++iteration) {
				long alignmentIndex = minAlignmentLongRead.size() - std::min(iteration + 1, (long)minAlignmentLongRead.size());
				pairedDBG.loadLongReadLink(minAlignmentLongRead[alignmentIndex], LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, this->contigMaxK, numThread);

				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());
					pairedDBG.setTolerence(std::min(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * longReadLibraryMT[0].getAverageInsSize(), 0.5 * pairedDBG.getCutoffLength()));
					pairedDBG.setMinLink(minLink);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
					pairedDBG.deleteErroneousEdgeScore(0.125, numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffoldCombine();
				}

				if (divisionFlag)
					pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
			}
		}

		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();


		if (libraryMT.size() > 0) {
			for (long iteration = 0; iteration < numIteration; ++iteration) {
				if (longReadLibraryMT.size() > 0) {
					long alignmentIndex = minAlignmentLongRead.size() - std::min(iteration + 1, (long)minAlignmentLongRead.size());
					pairedDBG.loadLongReadLink(minAlignmentLongRead[alignmentIndex], LONG_READ_MIN_ALIGNMENT_COVERAGE, minIdentity, this->contigMaxK, numThread);
				}

				pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
				long linkThreashold = (1 + iteration) * minLink;
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
				pairedDBG.setMinLink(linkThreashold);
				pairedDBG.setCutoffLength(2.0 * MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();

				if (numIteration >= 4 && outerIteration < 2)
					pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);

				if (divisionFlag)
					pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
			}
		}


		pairedDBG.trimSparseEnd();
		pairedDBG.trimRepeatEnd();
	}

}

void SolveDBG::execMinialign(const string targetFilename, const string &readFilename, const string &outFilename, const long numThread, const string minialignExecutable)
{
	std::ostringstream oss;

	oss << minialignExecutable <<  " -x pacbio -m 0 -O paf" << " -t " << numThread << " " << targetFilename << " " << readFilename << " >" << outFilename;

	std::cerr << "Executing minialign ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minialign fineshed." << endl;
}

void SolveDBG::execMinimap(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long minAlignmentLength, const long numThread, const string minimapExecutable)
{
	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimapExecutable <<  " -t " << numThread <<  " -L " << minAlignmentLength << " - ";

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
		oss << " " << *itr;
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename;

	std::cerr << "Executing minimap ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minimap fineshed." << endl;
}

void SolveDBG::execEmem(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long minMatchLength, const long numThread, const string ememExecutable)
{
	std::ostringstream oss;

	oss << "export NUCMER_E_MEM_OUTPUT_DIRPATH=" << optionSingleArgs["-tmp"] << "; ";

	char contigTempFileName[platanus::globalTmpFileDir.size() + 7 + 1];
	strcpy(contigTempFileName, platanus::globalTmpFileDir.c_str());
	strcat(contigTempFileName, "/XXXXXX"); 
	int fd = mkstemp(contigTempFileName);
	if (fd == -1) {
		throw platanus::TMPError();
	}
	FILE *fp = fdopen(fd, "w+");
	fclose(fp);

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;
	oss << " >" << contigTempFileName << "; ";

	char readTempFileName[platanus::globalTmpFileDir.size() + 7 + 1];
	strcpy(readTempFileName, platanus::globalTmpFileDir.c_str());
	strcat(readTempFileName, "/XXXXXX"); 
	fd = mkstemp(readTempFileName);
	if (fd == -1) {
		throw platanus::TMPError();
	}
	fp = fdopen(fd, "w+");
	fclose(fp);

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
		platanus::FILECOMPRESSION format = platanus::checkFileCompression(*itr);
		if (format == platanus::FILECOMPRESSION::BZIP2)
			oss << "bzip2 -cd " << *itr << " >>" << readTempFileName << "; ";
		else if (format == platanus::FILECOMPRESSION::GZIP)
			oss << "gzip -cd " << *itr << " >>" << readTempFileName << "; ";
		else
			oss << "cat " << *itr << " >>" << readTempFileName << "; ";
	}

	oss << ememExecutable << " -b -n" << " -l " << minMatchLength << " -t " << numThread << " " << contigTempFileName << " " << readTempFileName;
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename << "; ";

	oss << "rm " << contigTempFileName << " " << readTempFileName;

	std::cerr << "Executing e-mem ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "e-mem finished." << endl;
}


void SolveDBG::execMinimap2(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long numThread, const string minimap2Executable, const string minimap2Option)
{
	if (readFilenames.empty())
		return;

	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimap2Executable << " -c " << " -t " << numThread;
	if (optionBool["-minimap2_sensitive"])
		oss << " -p 0";
	
	oss << " " << minimap2Option << " - ";


	bool bzip2Flag = false;
	char bzip2TempFileName[platanus::globalTmpFileDir.size() + 7 + 1];
	platanus::FILECOMPRESSION format;

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
		format = platanus::checkFileCompression(*itr);
		if (format == platanus::FILECOMPRESSION::BZIP2) {
			bzip2Flag = true;
			break;
		}
	}

	if (bzip2Flag) {
		strcpy(bzip2TempFileName, platanus::globalTmpFileDir.c_str());
		strcat(bzip2TempFileName, "/XXXXXX"); 

        int fd = mkstemp(bzip2TempFileName);
        if (fd == -1) {
            throw platanus::TMPError();
		}
        FILE *fp = fdopen(fd, "w+");
		fclose(fp);

		std::ostringstream bzip2Oss;
		bzip2Oss << "bzip2 -cd ";

		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
			format = platanus::checkFileCompression(*itr);
			if (format == platanus::FILECOMPRESSION::BZIP2)
				bzip2Oss << " " << *itr;
			else
				oss << " " << *itr;
		}

		bzip2Oss << " >"  << bzip2TempFileName;
		if (system(bzip2Oss.str().c_str()) != 0) {
			throw platanus::AlignerError();
		}

		oss << " " << bzip2TempFileName;
	}
	else {
		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
			oss << " " << *itr;
	}

	oss << " | perl -pne \'s/cg:Z:\\S+//\' ";
	oss << " >>" << outFilename;

	std::cerr << "Executing minimap2 ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	if (bzip2Flag) {
        unlink(bzip2TempFileName);
	}

	cerr << "minimap2 finished." << endl;
}


void SolveDBG::extendConsensusToEstimateInsertSize(void)
{
	pairedDBG.clearContigPreviousParentNodeID();

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);

	pairedDBG.makeGraph(numThread);
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();
	pairedDBG.joinUnambiguousNodePairIterative(numThread);


	pairedDBG.setMinOverlap(minOverlapForScaffolding);

	pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

	for (long iteration = 0; iteration < 2; ++iteration) {
		for (unsigned libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
			pairedDBG.setTargetLibraryIndex(libraryIndex);
			updateAndWriteInsertSize(libraryIndex);

			if (iteration == 0)
				pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
			else
				pairedDBG.setMinLink(minLink);

			cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.makeGraph(numThread);
				pairedDBG.setOppositeBubbleContigIDGapped(numThread);
				pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);

				pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
				if (iteration > 0)
					pairedDBG.deleteErroneousEdgeIterative(numThread);

				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();
			}
		}

		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::MIS, numThread);
	}

	pairedDBG.setMinOverlap(this->contigMaxK - 1);
}


void SolveDBG::generateGraphFastg()
{
	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.makeGraph(numThread);

	std::ostringstream outStream;
	outStream << optionSingleArgs["-o"] << "_overlap.fastg";
	pairedDBG.outputFastg(outStream.str());
	pairedDBG.clearEdges();


	for (unsigned i = 0; i < libraryMT.size(); ++i) {
		pairedDBG.setTargetLibraryIndex(i);

		pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE);
		if (i > 0 && libraryMT[i][0].getAverageInsSize() == 0) {
			vector<long> insSizeDistribution;
			libraryMT[i][0].readInsertSizeFile(insSizeDistribution);
			pairedDBG.insertSizeDistribution(libraryMT[i], insSizeDistribution, numThread);
			libraryMT[i][0].estimateInsSize(insSizeDistribution, libraryMT[i - 1][0].getAverageInsSize(), platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

			std::ostringstream outStream;
			outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_insFreq.tsv";
			printInsertSizeFreq(outStream.str(), insSizeDistribution);
		}
		cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

		unsigned tolerenceFactor = MAX_TOL_FACTOR;
		pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
		cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);

		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_paired_end.fastg";
		pairedDBG.outputFastg(outStream.str());
		}
		pairedDBG.clearEdges();


		pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE);
		pairedDBG.makeGraph(numThread);
		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_paired_end_overlap.fastg";
		pairedDBG.outputFastg(outStream.str());
		}
		pairedDBG.clearEdges();
	}


	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
		pairedDBG.setTolerence(2 * this->contigMaxK);
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);
		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_long_read.fastg";
		pairedDBG.outputFastg(outStream.str());
		pairedDBG.clearEdges();
		}

		pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
		pairedDBG.makeGraph(numThread);
		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_long_read_overlap.fastg";
		pairedDBG.outputFastg(outStream.str());
		}
		pairedDBG.clearEdges();
	}
}


void SolveDBG::readLongReadPAFfileAndSaveLink(const string PAFFilename, platanus::Contig &contig, vector<SeqLib> &library, const long minAlignmentLength, const double minCoverage, const double minIdentity, const long tolerenceLength, const long numThread)
{
    cerr << "reading result of aligner ..." << endl;

	std::string oneLine, preQName(""), qName, tName, strand;
	long qLength, qStart, qEnd, tLength, tStart, tEnd, match, readLength=0;
	long threadIndex = 0, sumQAlignment = 0, sumTAlignment = 0;

	std::ifstream ifs;

	ifs.open(PAFFilename);
    while (ifs && getline(ifs, oneLine)) {
		std::istringstream ss(oneLine);
		ss >> qName >> qLength >> qStart >> qEnd >> strand >> tName >> tLength >> tStart >> tEnd >> match;
		sumQAlignment += qEnd - qStart;
		sumTAlignment += tEnd - tStart;
	}
	ifs.close();

	double readInsertionRate = (double)sumQAlignment / sumTAlignment;

	vector<platanus::LongReadAlignment> alignmentBuffer;
    long numMapped = 0;
    long sumMappedReadLength = 0;

    for (unsigned i = 0; i < numThread; ++i) {
		if (library[i].alignmentFP != NULL)
			fclose(library[i].alignmentFP);
		library[i].alignmentFP = platanus::makeTemporaryFile();
	}

	if (library[0].insertLengthFP != NULL)
		fclose(library[0].insertLengthFP);
	library[0].insertLengthFP = platanus::makeTemporaryFile();


	ifs.open(PAFFilename);
	while (1) {
		getline(ifs, oneLine);
		std::istringstream ss(oneLine);
		ss >> qName >> qLength >> qStart >> qEnd >> strand >> tName >> tLength >> tStart >> tEnd >> match;

		std::string tag;
		long score = 0;
		while (!ss.eof()) {
			ss >> tag;
			if (tag.size() > 5 && tag.substr(0, 5) == "AS:i:") {
				score = std::stoi(tag.substr(5));
				break;
			}
		}

		if (ifs.eof() || preQName != qName) {
			if (!alignmentBuffer.empty()) {
				fwrite(&readLength, sizeof(long), 1, library[threadIndex].alignmentFP);

				size_t bufferSize = alignmentBuffer.size();
				fwrite(&bufferSize, sizeof(size_t), 1, library[threadIndex].alignmentFP);
				for (unsigned j = 0; j < alignmentBuffer.size(); ++j)
					fwrite(&(alignmentBuffer[j]), sizeof(platanus::LongReadAlignment), 1, library[threadIndex].alignmentFP);

				threadIndex = (threadIndex + 1) %  numThread;

				fwrite(&readLength, sizeof(long), 1, library[0].insertLengthFP);
			}

			alignmentBuffer.clear();
			alignmentBuffer.clear();
			++numMapped;
			sumMappedReadLength += readLength;

			if (ifs.eof())
				break;
		}

		int alignmentLength = std::max(qEnd - qStart, tEnd - tStart);

		if ((double)match / alignmentLength >= minIdentity && (alignmentLength >= minAlignmentLength || (double)alignmentLength / std::min(qLength, tLength) >= minCoverage)) {
			unsigned contigIndex = contig.nameIndex[tName];
			alignmentBuffer.resize(alignmentBuffer.size() + 1);

			alignmentBuffer.back().score = score;
			alignmentBuffer.back().match = match;
			alignmentBuffer.back().readStart = qStart;
			alignmentBuffer.back().readEnd = qEnd;
			if (strand == "+") {
				alignmentBuffer.back().targetPosition.id = contigIndex + 1;
				alignmentBuffer.back().targetPosition.offset = tStart - qStart/readInsertionRate;
				alignmentBuffer.back().targetStart = tStart;
				alignmentBuffer.back().targetEnd = tEnd;
			}
			else {
				alignmentBuffer.back().targetPosition.id = -(contigIndex + 1);
				alignmentBuffer.back().targetPosition.offset = (tEnd - 1) + qStart/readInsertionRate;
				alignmentBuffer.back().targetStart = tEnd - 1;
				alignmentBuffer.back().targetEnd = tStart - 1;
			}

		}

		readLength = qLength;
		preQName = qName;
	}
	ifs.close();

	long contigSeqSize = 0;
    for (auto it = contig.seq.begin(); it != contig.seq.end(); ++it)
       contigSeqSize += it->base.size();

	library[0].setAverageCoverage(static_cast<double>(sumMappedReadLength) / contigSeqSize);

    cerr << "MAPPED_READ = " << numMapped << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}
