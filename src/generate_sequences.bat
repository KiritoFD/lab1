@echo off
echo DNA Sequence Generator with Answer Key
echo ====================================
echo.

set /p ref_length=Enter reference sequence length (default 100000): 
set /p query_length=Enter query sequence length (default 100000): 

if "%ref_length%"=="" set ref_length=100000
if "%query_length%"=="" set query_length=100000

echo.
echo Generating sequences with:
echo - Reference length: %ref_length% bp
echo - Query length: %query_length% bp
echo.

python generate_sequences.py -ref_length %ref_length% -query_length %query_length%

echo.
echo Additional options:
echo python generate_sequences.py -ref_length 50000 -query_length 100000 -ref custom_ref.txt -query custom_query.txt -answers custom_answers.txt
pause
