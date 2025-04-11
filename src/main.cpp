#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "aligner.h"
#include "fasta_reader.h"
#include "params.h"
#include "logger.h"

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " [OPTIONS] <reference.fa> <query.fa>\n"
              << "Options:\n"
              << "  -o, --output FILE       Output file (default: stdout)\n"
              << "  -m, --match INT         Match score (default: 2)\n"
              << "  -x, --mismatch INT      Mismatch penalty (default: -3)\n"
              << "  -g, --gap-open INT      Gap open penalty (default: -5)\n"
              << "  -e, --gap-extend INT    Gap extension penalty (default: -1)\n"
              << "  -l, --min-length INT    Minimum segment length (default: 50)\n"
              << "  -s, --min-score INT     Minimum alignment score (default: 30)\n"
              << "  -v, --verbose           Verbose output\n"
              << "  -h, --help              Print this help message\n";
}

int main(int argc, char* argv[]) {
    Logger logger;
    AlignmentParams params;
    std::string refFile, queryFile, outputFile;
    bool outputToFile = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-v" || arg == "--verbose") {
            logger.setVerbosity(Logger::INFO);
        } else if (arg == "-m" || arg == "--match") {
            if (i + 1 < argc) {
                params.matchScore = std::stoi(argv[++i]);
            }
        } else if (arg == "-x" || arg == "--mismatch") {
            if (i + 1 < argc) {
                params.mismatchPenalty = std::stoi(argv[++i]);
            }
        } else if (arg == "-g" || arg == "--gap-open") {
            if (i + 1 < argc) {
                params.gapOpenPenalty = std::stoi(argv[++i]);
            }
        } else if (arg == "-e" || arg == "--gap-extend") {
            if (i + 1 < argc) {
                params.gapExtendPenalty = std::stoi(argv[++i]);
            }
        } else if (arg == "-l" || arg == "--min-length") {
            if (i + 1 < argc) {
                params.minSegmentLength = std::stoi(argv[++i]);
            }
        } else if (arg == "-s" || arg == "--min-score") {
            if (i + 1 < argc) {
                params.minAlignmentScore = std::stoi(argv[++i]);
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                outputFile = argv[++i];
                outputToFile = true;
            }
        } else if (arg[0] != '-' && refFile.empty()) {
            refFile = arg;
        } else if (arg[0] != '-' && queryFile.empty()) {
            queryFile = arg;
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    if (refFile.empty() || queryFile.empty()) {
        std::cerr << "Error: Reference and query files must be specified.\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Set up output stream
    std::ofstream outFile;
    std::ostream* out = &std::cout;
    if (outputToFile) {
        outFile.open(outputFile);
        if (!outFile) {
            std::cerr << "Error: Cannot open output file " << outputFile << std::endl;
            return 1;
        }
        out = &outFile;
    }
    
    logger.log(Logger::INFO, "Starting SV-Aware Sequence Aligner");
    logger.log(Logger::INFO, "Reference file: " + refFile);
    logger.log(Logger::INFO, "Query file: " + queryFile);
    
    try {
        // Read reference sequences
        FastaReader refReader(refFile);
        std::vector<Sequence> references = refReader.readAll();
        logger.log(Logger::INFO, "Loaded " + std::to_string(references.size()) + " reference sequences");
        
        // Create and initialize aligner
        Aligner aligner(references, params, logger);
        
        // Read and align query sequences
        FastaReader queryReader(queryFile);
        std::vector<Sequence> queries = queryReader.readAll();
        logger.log(Logger::INFO, "Processing " + std::to_string(queries.size()) + " query sequences");
        
        // Write TSV header
        *out << "q_name\tq_len\tq_st\tq_en\tr_name\tr_len\tr_st\tr_en\tstrand\tscore\tedit_distance\tcigar\n";
        
        // Process each query
        for (const auto& query : queries) {
            logger.log(Logger::INFO, "Aligning query: " + query.name);
            std::vector<AlignmentResult> results = aligner.align(query);
            
            // Output results in TSV format
            for (const auto& result : results) {
                *out << result.toString() << "\n";
            }
        }
        
        logger.log(Logger::INFO, "Alignment completed successfully");
        
    } catch (const std::exception& e) {
        logger.log(Logger::ERROR, std::string("Error: ") + e.what());
        return 1;
    }
    
    if (outputToFile) {
        outFile.close();
    }
    
    return 0;
}
