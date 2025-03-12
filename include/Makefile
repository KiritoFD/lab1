# Makefile for DNA Repeat Finder project

# Compiler settings
CC = gcc
CFLAGS = -Wall -Wextra -fopenmp -O2
LDFLAGS = -fopenmp -lm

# Debug configuration
DEBUG_CFLAGS = -g -DDEBUG

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
INCLUDE_DIR = include

# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SOURCES))
DEPS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.d, $(SOURCES))

# Target executable
TARGET = $(BIN_DIR)/dna_repeat_finder

# Make all
all: directories $(TARGET)

# Debug build
debug: CFLAGS += $(DEBUG_CFLAGS)
debug: all

# Create necessary directories
directories:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR) $(INCLUDE_DIR)

# Link object files to create the executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -MMD -c $< -o $@

# Include dependency files
-include $(DEPS)

# Clean up
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.d $(TARGET)

# Clean everything including the compiled program
distclean: clean
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Install the program (to /usr/local/bin by default)
install: all
	install -m 755 $(TARGET) /usr/local/bin/

# Phony targets
.PHONY: all debug clean distclean directories install

# Help message
help:
	@echo "DNA Repeat Finder Makefile"
	@echo "Available targets:"
	@echo "  all       - Build the DNA repeat finder (default)"
	@echo "  debug     - Build with debug information"
	@echo "  clean     - Remove object files"
	@echo "  distclean - Remove object files and binary"
	@echo "  install   - Install to /usr/local/bin"
	@echo "  help      - Show this help message"
