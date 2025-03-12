# Makefile for DNA Repeat Finder project - Optimized for WSL/Linux

# Check if running on WSL for specific optimizations
WSL_CHECK := $(shell grep -q "Microsoft" /proc/version 2>/dev/null && echo 1 || echo 0)

# Compiler settings with specific optimizations
CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra -fopenmp -Iinclude -march=native -mtune=native \
         -O3 -ffast-math -funroll-loops -flto -fomit-frame-pointer -mavx2 -mfma -mbmi2 -mpopcnt \
         -fprefetch-loop-arrays -ftree-vectorize -fgraphite-identity -floop-interchange -floop-strip-mine -floop-block

# Use specific Zen4 optimizations if explicitly specified
ifeq ($(ARCH),zen4)
CFLAGS += -march=znver4 -mtune=znver4
endif

# WSL-specific optimizations to improve I/O performance
ifeq ($(WSL_CHECK),1)
CFLAGS += -DWSL_BUILD -pipe
endif

LDFLAGS = -fopenmp -lm -flto -fuse-linker-plugin

# Debug configuration
DEBUG_CFLAGS = -g -DDEBUG

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
INCLUDE_DIR = include
CPP_DIR = cpp

# Define sources and objects for main application
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
SRC_CORE_FILES := $(wildcard $(SRC_DIR)/core/*.c)
ALL_SRC_FILES := $(SRC_FILES) $(SRC_CORE_FILES)

# Exclude specific file from compilation
EXCLUDE_FILES := $(SRC_DIR)/dna_repeat_finder.c

# Filter out excluded files from source files
ALL_SRC_FILES := $(filter-out $(EXCLUDE_FILES), $(SRC_FILES) $(SRC_CORE_FILES))

# Main source file
MAIN_SOURCE = $(SRC_DIR)/main.c

# Exclude main source file from objects compilation
OBJ_FILES := $(filter-out $(OBJ_DIR)/main.o, $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES)))
OBJ_CORE_FILES := $(patsubst $(SRC_DIR)/core/%.c, $(OBJ_DIR)/core/%.o, $(SRC_CORE_FILES))
ALL_OBJ_FILES := $(OBJ_FILES) $(OBJ_CORE_FILES)

# Update object files list
OBJ_FILES := $(filter-out $(OBJ_DIR)/dna_repeat_finder.o, $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(ALL_SRC_FILES)))

# Main object file
MAIN_OBJECT = $(OBJ_DIR)/main.o

# Target executable
MAIN_TARGET = $(BIN_DIR)/dna_repeat_finder

# CPU core count for parallel compilation
CPU_CORES = $(shell nproc)

# Make all targets
all: directories $(MAIN_TARGET)

# Debug build
debug: CFLAGS += $(DEBUG_CFLAGS)
debug: all

# Build optimized for throughput
throughput: CFLAGS += -DOPTIMIZE_THROUGHPUT
throughput: all

# Create necessary directories
directories:
	@mkdir -p $(OBJ_DIR) $(OBJ_DIR)/core $(BIN_DIR)

# Link main application
$(MAIN_TARGET): $(ALL_OBJ_FILES) $(MAIN_OBJECT)
	$(CC) $^ -o $@ $(LDFLAGS)

# Compile regular source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Compile core source files
$(OBJ_DIR)/core/%.o: $(SRC_DIR)/core/%.c
	@mkdir -p $(OBJ_DIR)/core
	$(CC) $(CFLAGS) -MMD -c $< -o $@

# Dependencies
-include $(OBJ_DIR)/*.d
-include $(OBJ_DIR)/core/*.d

# Parallel compilation
.PHONY: parallel
parallel:
	@make -j$(CPU_CORES)

# Clean up
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d $(OBJ_DIR)/core/*.o $(OBJ_DIR)/core/*.d

# Clean everything including the compiled program
distclean: clean
	rm -f $(MAIN_TARGET)
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all debug throughput clean distclean directories parallel

# Default to parallel build for better performance
default: parallel

# Help message
help:
	@echo "DNA Repeat Finder Makefile - Optimized Build System"
	@echo "Available targets:"
	@echo "  all         - Build all applications"
	@echo "  parallel    - Build using all CPU cores ($(CPU_CORES) cores detected)"
	@echo "  debug       - Build with debug information"
	@echo "  clean       - Remove object files"
	@echo "  distclean   - Remove object files and binaries"
	@echo "  help        - Show this help message"
	@echo ""
	@echo "Special environment options:"
	@echo "  ARCH=zen4   - Enable specific AMD Zen4 optimizations"
