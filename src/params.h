#pragma once

#include <string>

/**
 * Structure to hold alignment parameters
 */
struct AlignmentParams {
    int matchScore = 2;
    int mismatchPenalty = -3;
    int gapOpenPenalty = -5;
    int gapExtendPenalty = -1;
    int minSegmentLength = 50;
    int minAlignmentScore = 30;
    int kmerSize = 15;  // Size of k-mers for seeding
    int minChainScore = 20;  // Minimum chain score to consider
};

/**
 * Structure to represent a DNA sequence
 */
class Sequence {
public:
    std::string name;
    std::string data;
    
    Sequence() = default;
    
    Sequence(const std::string& name, const std::string& data) 
        : name(name), data(data) {}
    
    size_t length() const {
        return data.length();
    }
    
    // Get reverse complement of the sequence
    Sequence reverseComplement() const {
        std::string revComp = data;
        for (size_t i = 0; i < revComp.size(); ++i) {
            char c = std::toupper(revComp[i]);
            switch (c) {
                case 'A': revComp[i] = 'T'; break;
                case 'C': revComp[i] = 'G'; break;
                case 'G': revComp[i] = 'C'; break;
                case 'T': revComp[i] = 'A'; break;
                default:  revComp[i] = 'N'; break;
            }
        }
        std::reverse(revComp.begin(), revComp.end());
        return Sequence(name, revComp);
    }
};

/**
 * Structure to represent a seed (match) between query and reference
 */
struct Seed {
    int queryPos;
    int refPos;
    int length;
    char strand;  // '+' or '-'
    
    Seed(int qp, int rp, int len, char str = '+') 
        : queryPos(qp), refPos(rp), length(len), strand(str) {}
    
    // For sorting seeds by query position
    bool operator<(const Seed& other) const {
        return queryPos < other.queryPos;
    }
};

/**
 * Structure to represent a chain of seeds
 */
struct Chain {
    std::vector<Seed> seeds;
    int score;
    int queryStart;
    int queryEnd;
    int refStart;
    int refEnd;
    std::string refName;
    char strand;
    
    Chain() : score(0), queryStart(0), queryEnd(0), refStart(0), refEnd(0), strand('+') {}
};

/**
 * Structure to represent an alignment result
 */
class AlignmentResult {
public:
    std::string queryName;
    int queryLength;
    int queryStart;
    int queryEnd;
    std::string refName;
    int refLength;
    int refStart;
    int refEnd;
    char strand;
    int score;
    int editDistance;
    std::string cigar;
    
    // Format result as TSV line
    std::string toString() const {
        return queryName + "\t" +
               std::to_string(queryLength) + "\t" +
               std::to_string(queryStart) + "\t" +
               std::to_string(queryEnd) + "\t" +
               refName + "\t" +
               std::to_string(refLength) + "\t" +
               std::to_string(refStart) + "\t" +
               std::to_string(refEnd) + "\t" +
               std::string(1, strand) + "\t" +
               std::to_string(score) + "\t" +
               std::to_string(editDistance) + "\t" +
               cigar;
    }
};
