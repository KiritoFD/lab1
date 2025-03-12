@echo off
echo DNA Repeat Finder
echo ================
echo.

if "%~2"=="" (
    echo Usage: %0 reference_file query_file
    echo Example: %0 reference.txt query.txt
    goto :EOF
)

echo Running DNA repeat finder with:
echo Reference file: %1
echo Query file: %2
echo.

dna_repeat_finder.exe %1 %2

echo.
echo Analysis complete.
pause
