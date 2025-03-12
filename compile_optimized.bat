@echo off
echo Compiling optimized DNA repeat finder with thread control...
gcc -O3 -fopenmp dna_repeat_finder.c -o dna_repeat_finder_optimized.exe
echo Compilation complete. Run dna_repeat_finder_optimized.exe to start.
pause
