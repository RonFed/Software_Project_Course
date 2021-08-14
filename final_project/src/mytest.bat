
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans
@REM .\spkmeans.exe 3 jacobi ..\tests_data\my_tests\data3_data.csv
@REM .\spkmeans.exe 3 ddg ..\..\data1_data.csv
@REM .\spkmeans.exe 3 jacobi ..\..\data1_data.csv
@REM .\spkmeans.exe 3 jacobi ..\..\data2_data.csv
@REM .\spkmeans.exe 0 spk ..\..\data2_data.csv
@REM .\spkmeans.exe 3 jacobi rami.txt
.\spkmeans.exe 0 wam ..\tests_data\elad_tests_general\input_5.txt
@REM .\spkmeans.exe 0 jacobi ..\tests_data\elad_tests_jacobi_2\input_J_4.txt