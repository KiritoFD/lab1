CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -I./src
LDFLAGS = 

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# Source files
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# Main target
TARGET = $(BIN_DIR)/sv_aligner

.PHONY: all clean directories

all: directories $(TARGET)

directories:
	mkdir -p $(BUILD_DIR) $(BIN_DIR)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
