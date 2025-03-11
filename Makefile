CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

all: dna_repeat_finder

dna_repeat_finder: dna_repeat_finder.cpp
	$(CXX) $(CXXFLAGS) -o dna_repeat_finder dna_repeat_finder.cpp

clean:
	rm -f dna_repeat_finder
