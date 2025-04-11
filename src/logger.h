#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

/**
 * Simple logger class for the aligner
 */
class Logger {
public:
    enum LogLevel { 
        ERROR = 0,    // Fatal errors
        WARNING = 1,  // Non-fatal issues
        INFO = 2,     // Main program steps
        DEBUG = 3     // Detailed debug information
    };

private:
    LogLevel verbosity;

public:
    Logger() : verbosity(WARNING) {}
    
    void setVerbosity(LogLevel level) {
        verbosity = level;
    }
    
    void log(LogLevel level, const std::string& message) {
        if (level <= verbosity) {
            auto now = std::chrono::system_clock::now();
            auto now_time_t = std::chrono::system_clock::to_time_t(now);
            std::tm tm = *std::localtime(&now_time_t);
            
            std::cerr << "[" << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "] ";
            
            switch (level) {
                case ERROR:   std::cerr << "ERROR: "; break;
                case WARNING: std::cerr << "WARNING: "; break;
                case INFO:    std::cerr << "INFO: "; break;
                case DEBUG:   std::cerr << "DEBUG: "; break;
            }
            
            std::cerr << message << std::endl;
        }
    }
};
