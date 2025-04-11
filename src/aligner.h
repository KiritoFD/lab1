#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <limits>
#include <cmath>
#include "params.h"
#include "logger.h"

/**
 * Index for efficient k-mer lookup in reference sequences
 */
class ReferenceIndex {
private:
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> kmerIndex; // kmer -> [(refId, pos)]
    std::vector<Sequence> references;
    int kmerSize;
    
public:
    ReferenceIndex(const std::vector<Sequence>& refs, int kSize) 
        : references(refs), kmerSize(kSize) {
        buildIndex();
    }
    
    void buildIndex() {
        for (int refId = 0; refId < references.size(); ++refId) {
            const auto& ref = references[refId];
            for (int i = 0; i <= ref.data.length() - kmerSize; ++i) {
                std::string kmer = ref.data.substr(i, kmerSize);
                kmerIndex[kmer].emplace_back(refId, i);
            }
        }
    }
    
    // Find all occurrences of a k-mer in the reference
    std::vector<std::pair<int, int>> lookup(const std::string& kmer) const {
        auto it = kmerIndex.find(kmer);
        if (it != kmerIndex.end()) {
            return it->second;
        }
        return {};
    }
    
    const std::vector<Sequence>& getReferences() const {
        return references;
    }
};

/**
 * Class implementing the SV-aware sequence alignment algorithm
 */
class Aligner {
private:
    ReferenceIndex refIndex;
    AlignmentParams params;
    Logger& logger;
    
public:
    Aligner(const std::vector<Sequence>& references, const AlignmentParams& p, Logger& log) 
        : refIndex(references, p.kmerSize), params(p), logger(log) {}
    
    std::vector<AlignmentResult> align(const Sequence& query) {
        logger.log(Logger::DEBUG, "Finding seeds for query: " + query.name);
        
        // Find seeds (exact k-mer matches) between query and references
        std::vector<Seed> forwardSeeds = findSeeds(query, '+');
        std::vector<Seed> reverseSeeds = findSeeds(query.reverseComplement(), '-');
        
        // Combine and sort all seeds by query position
        std::vector<Seed> allSeeds = forwardSeeds;
        allSeeds.insert(allSeeds.end(), reverseSeeds.begin(), reverseSeeds.end());
        std::sort(allSeeds.begin(), allSeeds.end());
        
        logger.log(Logger::DEBUG, "Found " + std::to_string(allSeeds.size()) + " seeds");
        
        // Chain seeds into consistent alignments
        std::vector<Chain> chains = chainSeeds(allSeeds, query);
        logger.log(Logger::DEBUG, "Created " + std::to_string(chains.size()) + " seed chains");
        
        // Refine chains to get detailed alignments
        std::vector<AlignmentResult> results;
        for (const auto& chain : chains) {
            if (chain.score >= params.minChainScore) {
                AlignmentResult result = refineAlignment(chain, query);
                if (result.score >= params.minAlignmentScore && 
                    result.queryEnd - result.queryStart >= params.minSegmentLength) {
                    results.push_back(result);
                }
            }
        }
        
        // Remove overlapping segments in query sequence
        results = removeOverlappingResults(results);
        
        logger.log(Logger::DEBUG, "Final alignment count: " + std::to_string(results.size()));
        return results;
    }

private:
    // Find exact k-mer matches between query and reference
    std::vector<Seed> findSeeds(const Sequence& query, char strand) {
        std::vector<Seed> seeds;
        const int kSize = params.kmerSize;
        
        for (int i = 0; i <= query.length() - kSize; ++i) {
            std::string kmer = query.data.substr(i, kSize);
            auto matches = refIndex.lookup(kmer);
            
            for (const auto& match : matches) {
                int refId = match.first;
                int refPos = match.second;
                seeds.emplace_back(i, refPos, kSize, strand);
                seeds.back().refName = refIndex.getReferences()[refId].name;
            }
        }
        
        return seeds;
    }
    
    // Chain seeds into consistent alignments
    std::vector<Chain> chainSeeds(const std::vector<Seed>& seeds, const Sequence& query) {
        std::vector<Chain> chains;
        // Implementation of a seed chaining algorithm
        // This is a simplified version that would need to be expanded
        
        if (seeds.empty()) return chains;
        
        // Group seeds by reference sequence and strand
        std::unordered_map<std::string, std::vector<Seed>> seedsByRef;
        for (const auto& seed : seeds) {
            std::string key = seed.refName + seed.strand;
            seedsByRef[key].push_back(seed);
        }
        
        // Chain seeds for each reference/strand combination
        for (const auto& entry : seedsByRef) {
            const auto& refSeeds = entry.second;
            
            // Simple greedy chaining as a placeholder
            // A more sophisticated dynamic programming approach would be used in practice
            Chain currentChain;
            currentChain.refName = refSeeds[0].refName;
            currentChain.strand = refSeeds[0].strand;
            
            for (const auto& seed : refSeeds) {
                // Add to current chain if it extends consistently
                if (currentChain.seeds.empty() || 
                    (seed.queryPos >= currentChain.seeds.back().queryPos + currentChain.seeds.back().length &&
                     seed.refPos >= currentChain.seeds.back().refPos + currentChain.seeds.back().length)) {
                    
                    currentChain.seeds.push_back(seed);
                    currentChain.score += seed.length;
                    
                    // Update chain boundaries
                    if (currentChain.seeds.size() == 1) {
                        currentChain.queryStart = seed.queryPos;
                        currentChain.refStart = seed.refPos;
                    }
                    currentChain.queryEnd = seed.queryPos + seed.length;
                    currentChain.refEnd = seed.refPos + seed.length;
                } 
                // If doesn't extend current chain, start a new one if current one is good enough
                else if (currentChain.score >= params.minChainScore) {
                    chains.push_back(currentChain);
                    
                    // Start new chain with this seed
                    currentChain = Chain();
                    currentChain.seeds.push_back(seed);
                    currentChain.score = seed.length;
                    currentChain.queryStart = seed.queryPos;
                    currentChain.refStart = seed.refPos;
                    currentChain.queryEnd = seed.queryPos + seed.length;
                    currentChain.refEnd = seed.refPos + seed.length;
                    currentChain.refName = seed.refName;
                    currentChain.strand = seed.strand;
                }
            }
            
            // Add the final chain if it's good enough
            if (currentChain.score >= params.minChainScore) {
                chains.push_back(currentChain);
            }
        }
        
        return chains;
    }
    
    // Refine a chain to get detailed alignment
    AlignmentResult refineAlignment(const Chain& chain, const Sequence& query) {
        AlignmentResult result;
        result.queryName = query.name;
        result.queryLength = query.length();
        result.queryStart = chain.queryStart;
        result.queryEnd = chain.queryEnd;
        result.refName = chain.refName;
        
        // Find the reference sequence
        const Sequence* refSeq = nullptr;
        for (const auto& ref : refIndex.getReferences()) {
            if (ref.name == chain.refName) {
                refSeq = &ref;
                break;
            }
        }
        
        if (!refSeq) {
            throw std::runtime_error("Reference sequence not found: " + chain.refName);
        }
        
        result.refLength = refSeq->length();
        result.refStart = chain.refStart;
        result.refEnd = chain.refEnd;
        result.strand = chain.strand;
        
        // For now, use a simple scoring approach based on chain score
        // In a real implementation, we would run Smith-Waterman here
        result.score = chain.score;
        result.editDistance = (result.queryEnd - result.queryStart) / 10; // Placeholder
        result.cigar = std::to_string(result.queryEnd - result.queryStart) + "M"; // Placeholder
        
        return result;
    }
    
    // Remove overlapping query segments, keeping highest scoring ones
    std::vector<AlignmentResult> removeOverlappingResults(std::vector<AlignmentResult> results) {
        if (results.size() <= 1) return results;
        
        // Sort by query name and then by score (descending)
        std::sort(results.begin(), results.end(), [](const AlignmentResult& a, const AlignmentResult& b) {
            if (a.queryName != b.queryName) return a.queryName < b.queryName;
            return a.score > b.score; // Higher score first
        });
        
        std::vector<AlignmentResult> filteredResults;
        std::vector<bool> keep(results.size(), true);
        
        // Mark overlapping segments for removal
        for (int i = 0; i < results.size(); i++) {
            if (!keep[i]) continue;
            
            for (int j = i + 1; j < results.size(); j++) {
                if (results[i].queryName != results[j].queryName) break; // Different query
                
                // Check for overlap
                if (results[i].queryStart < results[j].queryEnd && results[j].queryStart < results[i].queryEnd) {
                    keep[j] = false; // Mark the lower-scoring one for removal
                }
            }
        }
        
        // Create the filtered list
        for (int i = 0; i < results.size(); i++) {
            if (keep[i]) {
                filteredResults.push_back(results[i]);
            }
        }
        
        return filteredResults;
    }
};
