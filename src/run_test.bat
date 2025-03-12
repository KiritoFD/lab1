@echo off
echo DNA Repeat Finder Test Framework
echo ===============================
echo.

echo Compiling DNA repeat finder...
gcc -O2 -static dna_repeat_finder.c -o dna_repeat_finder.exe
if %errorlevel% neq 0 (
    echo Compilation failed!
    pause
    exit /b 1
)

echo.
echo Starting test sequence...
python test_repeat_finder.py

echo.
echo Test complete!
pause
