python generate_input.py 10 > my_intput_rand.txt
python kmeans.py 3 600 < tests/input_1.txt > my_output1py.txt
python kmeans.py 7 < tests/input_2.txt > my_output2py.txt
python kmeans.py 15 300 < tests/input_3.txt > my_output3py.txt
python kmeans.py 10 < my_intput_rand.txt > my_output_randpy.txt

gcc -ansi -Wall -Wextra -Werror -pedantic-errors kmeans.c -lm -o ex1
@REM gcc kmeans.c -o ex1

.\ex1.exe 3 600 < tests\input_1.txt > my_output1c.txt
.\ex1.exe 7 < tests\input_2.txt > my_output2c.txt
.\ex1.exe 15 300 < tests\input_3.txt > my_output3c.txt
.\ex1.exe 10 < my_intput_rand.txt > my_output_randc.txt

FC my_output1py.txt tests/output_1.txt
FC my_output2py.txt tests/output_2.txt
FC my_output3py.txt tests/output_3.txt

FC my_output1c.txt tests/output_1.txt
FC my_output2c.txt tests/output_2.txt
FC my_output3c.txt tests/output_3.txt

FC my_output_randpy.txt my_output_randc.txt
pause