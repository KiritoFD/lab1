#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include "params.h"

/**
 * Class implementing the Smith-Waterman local alignment algorithm with affine gap penalties
 */
class SmithWaterman {
private:
    const AlignmentParams& params;
    
public:
    SmithWaterman(const AlignmentParams& p) : params(p) {}
    
    // Structure to hold detailed alignment result
    struct SWResult {
        int score;
        int queryStart;
        int queryEnd;
        int refStart;
        int refEnd;
        int editDistance;
        std::string cigar;
    };
    
    // Perform local alignment between two sequences
    SWResult align(const std::string& query, const std::string& ref) {
        int n = query.length();
        int m = ref.length();
        
        // Initialize scoring matrices
        std::vector<std::vector<int>> score(n + 1, std::vector<int>(m + 1, 0));
        std::vector<std::vector<int>> E(n + 1, std::vector<int>(m + 1, 0));  // Gap in query
        std::vector<std::vector<int>> F(n + 1, std::vector<int>(m + 1, 0));  // Gap in reference
        std::vector<std::vector<char>> backtrack(n + 1, std::vector<char>(m + 1, 0));
        
        int maxScore = 0;
        int endI = 0, endJ = 0;
        
        // Fill in the DP matrices
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                // Compute E (gap in query)
                E[i][j] = std::max(E[i][j-1] + params.gapExtendPenalty, 
                                  score[i][j-1] + params.gapOpenPenalty + params.gapExtendPenalty);
                
                // Compute F (gap in reference)
                F[i][j] = std::max(F[i-1][j] + params.gapExtendPenalty,
                                  score[i-1][j] + params.gapOpenPenalty + params.gapExtendPenalty);
                
                // Match/mismatch score
                int match = score[i-1][j-1] + (query[i-1] == ref[j-1] ? params.matchScore : params.mismatchPenalty);
                
                // Compute the best score
                int bestScore = std::max({0, match, E[i][j], F[i][j]});
                score[i][j] = bestScore;
                
                // Record backtrace pointer
                if (bestScore == 0) backtrack[i][j] = 'S'; // Stop
                else if (bestScore == match) backtrack[i][j] = 'D'; // Diagonal
                else if (bestScore == E[i][j]) backtrack[i][j] = 'L'; // Left (gap in query)
                else backtrack[i][j] = 'U'; // Up (gap in reference)
                
                // Update global max
                if (bestScore > maxScore) {
                    maxScore = bestScore;
                    endI = i;
                    endJ = j;
                }
            }
        }
        
        // Trace back to find alignment
        SWResult result;
        result.score = maxScore;
        result.queryEnd = endI;
        result.refEnd = endJ;
        
        if (maxScore == 0) {
            // No alignment found
            result.queryStart = 0;
            result.refStart = 0;
            result.editDistance = 0;
            result.cigar = "0M";
            return result;
        }
        
        // Trace back to find the start of alignment
        int i = endI, j = endJ;
        int editDist = 0;
        std::string cigarOps;
        char currentOp = 'M';  // Start with match
        int opCount = 0;
        
        while (backtrack[i][j] != 'S') {
            char bt = backtrack[i][j];
            char newOp;
            
            if (bt == 'D') {
                // Diagonal move (match/mismatch)
                i--; j--;
                if (query[i] != ref[j]) editDist++;
                newOp = 'M';
            } else if (bt == 'L') {
                // Left move (gap in query)
                j--;
                editDist++;
                newOp = 'D';  // Deletion in the query
            } else { // 'U'
                // Up move (gap in reference)
                i--;
                editDist++;
                newOp = 'I';  // Insertion in the query
            }
            
            if (newOp != currentOp) {
                if (opCount > 0) {
                    cigarOps = std::to_string(opCount) + currentOp + cigarOps;
                }
                currentOp = newOp;
                opCount = 1;
            } else {
                opCount++;
            }
        }
        
        // Add the final CIGAR operation
        if (opCount > 0) {
            cigarOps = std::to_string(opCount) + currentOp + cigarOps;
        }
        
        result.queryStart = i;
        result.refStart = j;
        result.editDistance = editDist;
        result.cigar = cigarOps;
        
        return result;
    }
};
