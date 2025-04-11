#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include "params.h"

/**
 * Class for reading FASTA format sequence files
 */
class FastaReader {
private:
    std::string filename;

public:
    FastaReader(const std::string& filename) : filename(filename) {}
    
    // Read all sequences from a FASTA file
    std::vector<Sequence> readAll() {
        std::vector<Sequence> sequences;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        std::string line, name, data;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Process previous sequence if any
                if (!name.empty()) {
                    sequences.emplace_back(name, data);
                }
                
                // Start new sequence
                name = line.substr(1);
                // Remove any whitespace after the first space
                size_t spacePos = name.find_first_of(" \t");
                if (spacePos != std::string::npos) {
                    name = name.substr(0, spacePos);
                }
                data.clear();
            } else {
                // Normalize and append sequence data (convert to uppercase)
                for (char c : line) {
                    if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
                        data.push_back(std::toupper(c));
                    }
                }
            }
        }
        
        // Add the last sequence if any
        if (!name.empty()) {
            sequences.emplace_back(name, data);
        }
        
        return sequences;
    }
};
