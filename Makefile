# Makefile for DNA Repeat Finder project

# Compiler settings with extreme optimizations for R9 7940HX (Zen 4 architecture)
CC = gcc
CFLAGS = -Wall -Wextra -fopenmp -Iinclude \
         -march=znver4 -mtune=znver4 \
         -O3 -ffast-math -funroll-loops -flto \
         -fomit-frame-pointer -mfma -mavx2 -mbmi2 -mpopcnt \
         -fprefetch-loop-arrays -ftree-vectorize \
         -fgraphite-identity -floop-interchange -floop-strip-mine -floop-block

LDFLAGS = -fopenmp -lm -flto -fuse-linker-plugin

# Debug configuration
DEBUG_CFLAGS = -g -DDEBUG

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
INCLUDE_DIR = include

# Define all source files except for files with main() functions
COMMON_SOURCES = $(SRC_DIR)/dna_common.c \
                 $(SRC_DIR)/dna_io.c \
                 $(SRC_DIR)/dna_traditional.c \
                 $(SRC_DIR)/dna_graph.c \
                 $(SRC_DIR)/dna_sequence_utils.c

# Main application source
MAIN_SOURCE = $(SRC_DIR)/main.c

# Generator application source
GENERATOR_SOURCE = $(SRC_DIR)/generate_test_sequences.c

# Derive object files from sources
COMMON_OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(COMMON_SOURCES))
MAIN_OBJECT = $(OBJ_DIR)/main.o
GENERATOR_OBJECT = $(OBJ_DIR)/generate_test_sequences.o

# Target executables
MAIN_TARGET = $(BIN_DIR)/dna_repeat_finder
GENERATOR_TARGET = $(BIN_DIR)/generate_test_sequences

# CPU core count for parallel compilation
CPU_CORES = $(shell nproc)

# Make all
all: directories $(MAIN_TARGET) $(GENERATOR_TARGET)

# Debug build
debug: CFLAGS += $(DEBUG_CFLAGS)
debug: all

# Performance build (default)
performance: CFLAGS += -DNDEBUG
performance: all

# Create necessary directories
directories:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR) $(INCLUDE_DIR)

# Link main application
$(MAIN_TARGET): $(COMMON_OBJECTS) $(MAIN_OBJECT)
	$(CC) $^ -o $@ $(LDFLAGS)

# Link generator application
$(GENERATOR_TARGET): $(GENERATOR_OBJECT)
	$(CC) $^ -o $@ $(LDFLAGS)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Dependencies
-include $(OBJ_DIR)/*.d

# Parallel compilation using all available cores
.PHONY: parallel
parallel:
	@make -j$(CPU_CORES)

# Clean up
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d $(MAIN_TARGET) $(GENERATOR_TARGET)

# Clean everything including the compiled program
distclean: clean
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Install the program (to /usr/local/bin by default)
install: all
	install -m 755 $(MAIN_TARGET) /usr/local/bin/
	install -m 755 $(GENERATOR_TARGET) /usr/local/bin/

# Phony targets
.PHONY: all debug performance clean distclean directories install parallel

# Default to parallel build for better performance
default: parallel

# Help message
help:
	@echo "DNA Repeat Finder Makefile"
	@echo "Available targets:"
	@echo "  all         - Build all applications"
	@echo "  performance - Build with extreme R9 7940HX optimizations (default)"
	@echo "  parallel    - Build using all CPU cores ($(CPU_CORES) cores detected)"
	@echo "  debug       - Build with debug information"
	@echo "  clean       - Remove object files"
	@echo "  distclean   - Remove object files and binary"
	@echo "  install     - Install to /usr/local/bin"
	@echo "  help        - Show this help message"
