python ex1.py 3 600 < tests/input_1.txt > my_output1.txt
python ex1.py 7 < tests/input_2.txt > my_output2.txt
python ex1.py 15 300 < tests/input_3.txt > my_output3.txt

FC my_output1.txt tests/output_1.txt
FC my_output2.txt tests/output_2.txt
FC my_output3.txt tests/output_3.txt
pause