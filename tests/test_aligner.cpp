#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include "../src/aligner.h"
#include "../src/fasta_reader.h"
#include "../src/params.h"
#include "../src/logger.h"

// Test utilities
void assertEquals(const std::string& test_name, int expected, int actual) {
    if (expected != actual) {
        std::cerr << "FAIL: " << test_name << " - Expected " << expected << " but got " << actual << std::endl;
        exit(1);
    } else {
        std::cout << "PASS: " << test_name << std::endl;
    }
}

void assertEquals(const std::string& test_name, const std::string& expected, const std::string& actual) {
    if (expected != actual) {
        std::cerr << "FAIL: " << test_name << " - Expected \"" << expected << "\" but got \"" << actual << "\"" << std::endl;
        exit(1);
    } else {
        std::cout << "PASS: " << test_name << std::endl;
    }
}

// Test the basic alignment functionality
void testBasicAlignment() {
    Logger logger;
    logger.setVerbosity(Logger::DEBUG);
    
    // Create a simple reference sequence
    Sequence ref("ref1", "ACGTACGTACGTACGTACGT");
    std::vector<Sequence> references = {ref};
    
    // Create a query that matches exactly
    Sequence query("query1", "ACGTACGT");
    
    // Set up aligner with default parameters
    AlignmentParams params;
    Aligner aligner(references, params, logger);
    
    // Get alignment results
    std::vector<AlignmentResult> results = aligner.align(query);
    
    // Verify we got a result
    assertEquals("Basic alignment result count", 1, results.size());
    
    if (results.size() > 0) {
        // Check some basic properties of the alignment
        assertEquals("Query name", "query1", results[0].queryName);
        assertEquals("Reference name", "ref1", results[0].refName);
        assertEquals("Query start", 0, results[0].queryStart);
        assertEquals("Query end", 8, results[0].queryEnd);
        assertEquals("Strand", '+', results[0].strand);
    }
}

// Test alignment with inversions
void testInversionAlignment() {
    Logger logger;
    logger.setVerbosity(Logger::DEBUG);
    
    // Create reference sequence
    Sequence ref("ref1", "ACGTACGTACGTACGTACGT");
    std::vector<Sequence> references = {ref};
    
    // Create a query with a reverse complemented segment
    Sequence query("query1", "ACGTTACGT"); // ACGT (revcomp of ACGT) ACGT
    
    // Set up aligner with default parameters
    AlignmentParams params;
    params.minSegmentLength = 4; // Smaller for this test
    Aligner aligner(references, params, logger);
    
    // Get alignment results
    std::vector<AlignmentResult> results = aligner.align(query);
    
    // We expect two segments due to the inversion
    assertEquals("Inversion alignment result count", 2, results.size());
    
    if (results.size() >= 2) {
        // The first segment should be on + strand
        assertEquals("First segment strand", '+', results[0].strand);
        // The second segment should be on - strand (inversion)
        assertEquals("Second segment strand", '-', results[1].strand);
    }
}

// Add more tests as needed...

int main() {
    std::cout << "Running aligner tests...\n";
    
    testBasicAlignment();
    testInversionAlignment();
    
    std::cout << "All tests passed!\n";
    return 0;
}
