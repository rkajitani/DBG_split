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

#ifndef MAPPER_H
#define MAPPER_H

#include <climits>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include "common.h"
#include "seqlib.h"
#include "kmer.h"




//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// MapPointer class
// this class has kmer mapping position information for pool
// num means the number of kmer in pool exists
// MapPointer means the first position in pool
struct MapPointer
{
    unsigned num;
    platanus::Position *position;
    MapPointer(): num(0), position(NULL) {}
    ~MapPointer() {}

};
//////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////
// mapper class
//////////////////////////////////////////////////////////////////////////////////////
class Mapper
{
    typedef unsigned long long u64_t;

private:
	static const double MIN_IDENTITY_FOR_GAP_CLOSE;
	static const double MIN_IDENTITY_FOR_SCAFFOLD;
	static const double MIN_IDENTITY_TO_CHECK_MAPPING;
	static const int DEFAULT_MAPKEY_LIMIT;
	static const long DISPLAY_NUM_READ_UNIT;

    int keyLength;
	int mapKeyUpperLimit;
    u64_t mask;
    int seedLength;
	std::vector<int> multiSeedLength;
    int indexLength;
    long seqPoolSize;
    long numSeq;
    std::vector<platanus::SEQ> seq;
    std::unordered_map<u64_t, MapPointer> table;
    std::unique_ptr<platanus::Position[]> positionPool;

	std::unordered_map<std::string, unsigned> nameIndex;

    platanus::Position *reservePositionPool(const unsigned long long total);
    void fillPositionPool(FILE *tmpFP, platanus::Position *pos);
    void mapReadGapCloseLeft(const platanus::SEQ &read, const Kmer31 &kmer, std::vector<platanus::Position> &buffer);
    void mapReadGapCloseRight(const platanus::SEQ &read, const Kmer31 &kmer, std::vector<platanus::Position> &buffer);
    void mapReadGapClose(const platanus::SEQ &read, const Kmer31 &kmer, std::vector<platanus::Position> &buffer, const bool isLeft);
    long searchGapStart(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult);
    long searchGapEnd(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult, const long length);
    bool checkGapContinuious(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult, const long start, const long end, const long length);

    template <typename T>
    inline T sign(T val) { return (T(0) < val) - (T(0) > val); }


public:

    Mapper();
    Mapper(const long seed, const long key);
    Mapper(const Mapper &) = delete;
    Mapper &operator=(const Mapper &) = delete;
    ~Mapper() = default;



    // getter and setter
    int getSeedLength(void) const { return seedLength; }
    int getKeyLength(void) const { return keyLength; }
    long getSeqPoolSize(void) const { return seqPoolSize; }
    long getNumSeq(void) const { return numSeq; }
    const platanus::SEQ &getSeq(const long position) const { return seq.at(position); }
    long getSeqLength(const long position) const { return seq[position].length; }
    char getBase(const long seqID, const long position) const { return seq[seqID].base.at(position); }
    const std::vector<platanus::SEQ> &getSeqAll(void) const { return seq; }

    void setContig(const platanus::Contig &contig);
    void setKeyLength(const int key) {
		keyLength = key > 32 ? 32 : key;
		mask = keyLength < 32 ? ~(~0ull << (2 * keyLength)) : ~0ull;
	}
    void setSeedLength(const int inputSeedLength) { seedLength = inputSeedLength; }
	void setMultiSeedLength(const std::vector<int> &inputSeedLength) { multiSeedLength = inputSeedLength; }
	void setMapKeyUpperLimit(const int upperLimit) { mapKeyUpperLimit = upperLimit; }

    void releaseContig(platanus::Contig &contig) { contig.seq = std::move(seq); }
    void moveSeq(std::vector<platanus::SEQ> &sseq) { sseq = std::move(seq); }

    void insertHashKey(void);
    void makeKmerTable(void);
    long mapKey(const Kmer31 &kmer, std::vector<platanus::Position> &positionBuffer) const;
	long mapKeyWithLimit(const Kmer31 &kmer, std::vector<platanus::Position> &positionBuffer) const;
    platanus::Position mapRead(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer, const long wordLength);
	void mapReadMultiReports(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer, const long wordLength, const long minContigLength);
	void mapReadMultiReportsMultiSeed(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer, const long minContigLength);
	platanus::Position mapReadMultiSeedFiltered(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer);
	void mapReadMultiReportsMultiSeedFiltered(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer, const long minContigLength);
	platanus::Position mapReadUngapAlignment(const platanus::SEQ &read, const double minIdentity, std::vector<platanus::Position> &positionBuffer, double *returnIdentity);
	void mapReadUngapAlignmentMulti(const platanus::SEQ &read, const double minIdentity, std::vector<platanus::Position> &positionBuffer, std::vector<platanus::Position> &returnPositionBuffer);
	void filterMappedPositionUngapAlignment(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer, const double minIdentity);
	bool judgeMappedPositionUngapAlignment(const platanus::SEQ &read, const platanus::Position readLeftPosition, const double minIdentity);
    void mapPairMT(std::vector<SeqLib> &library, const long minInsertion, const long numThresod);
	void mapTagPairMT(std::vector<SeqLib> &library, const long numThread);
	void mapPairAndSaveReadLink(std::vector<SeqLib> &library, const long minInsertion, const long minContigLength, const long numThread);
    void mapSmallGap(std::vector<SeqLib> &library, FILE *gapFP, const long numThread);
    void gatherPairReadMappedSameContig(std::vector<SeqLib> &library, const long numThread);
	void mapPairToCalculateCoverage(std::vector<SeqLib> &library, const long numThread);
	void mapLongReadAndSaveReadLink(std::vector<SeqLib> &library, const long minContigLength, const long numThread);
	void mapAnchorBubbleMT(platanus::Contig &anchorContig, std::vector<platanus::Region> &anchorMap, const long numThread);
	void mapAnchorHomoMT(platanus::Contig &anchorContig, std::vector<platanus::Region> &anchorMap, const long numThread);
};



//////////////////////////////////////////////////////////////////////////////////////
// mapper class using bubble seq too
//////////////////////////////////////////////////////////////////////////////////////
class HeteroMapper
{
private:

	static const double MIN_IDENTITY_FOR_SCAFFOLD;

    long keyLength;
    long seedLength;
    std::vector<platanus::Position> bubblePosition;

public:

    Mapper contigMap;
    Mapper bubbleMap;

    HeteroMapper() = delete;
    HeteroMapper(const long seed, const long key): keyLength(key), seedLength(seed), bubblePosition(), contigMap(seed, key), bubbleMap(seed, key)
    { keyLength = std::min(key, static_cast<long>(32)); }
    HeteroMapper(const HeteroMapper &) = delete;
    HeteroMapper &operator=(const HeteroMapper &) = delete;
    ~HeteroMapper() = default;

    const std::vector<platanus::Position> &getBubblePositionRef(void) const
    {return bubblePosition; }

    void setContigMap(platanus::Contig &contig)
    {contigMap.setContig(contig); }
    void setBubbleMap(platanus::Contig &contig)
    {bubbleMap.setContig(contig); }
	void setMultiSeedLength(const std::vector<int> &inputSeedLength)
	{contigMap.setMultiSeedLength(inputSeedLength); bubbleMap.setMultiSeedLength(inputSeedLength); }
	void setKeyLength(const int key)
	{contigMap.setKeyLength(key); bubbleMap.setKeyLength(key); }

    void makeKmerTableContigMap(void)
    {contigMap.makeKmerTable(); }
    void makeKmerTableBubbleMap(void)
    {bubbleMap.makeKmerTable(); }

    void mergeBubble(const long numThread);
    void mapPairMT(std::vector<SeqLib> &library, const long minInsertion, const long numThread);
    platanus::Position mapRead(const platanus::SEQ &read, std::vector<platanus::Position> &positionBuffer, const long wordLength);

};



#endif
