@echo off
echo Compiling optimized DNA sequence generator for 100,000 bp sequences...
gcc -O3 generate_test_sequences.c -o generate_sequences.exe
echo Compilation complete.
echo.
echo Usage examples:
echo   generate_sequences.exe                    - Creates 100,000 bp sequences (default)
echo   generate_sequences.exe -length 200000     - Creates 200,000 bp sequences
echo   generate_sequences.exe -ref reference.txt -query query.txt  - Custom filenames
pause
