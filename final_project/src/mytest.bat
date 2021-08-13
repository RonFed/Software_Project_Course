
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans

.\spkmeans.exe 3 wam ..\tests_data\my_tests\data1_data.csv
@REM .\spkmeans.exe 3 ddg ..\..\data1_data.csv
@REM .\spkmeans.exe 3 jacobi ..\..\data1_data.csv
@REM .\spkmeans.exe 3 jacobi ..\..\data2_data.csv
@REM .\spkmeans.exe 0 spk ..\..\data2_data.csv
@REM .\spkmeans.exe 3 jacobi rami.txt
.\spkmeans.exe 0 lnorm ..\tests_data\elad_tests_general\input_1.txt