Compile:
 g++ main.cpp -o main.exe

Run (in cmd):
 main.exe < ($testFileName).txt

testFileName:
 brandy
 dualtest
 emptyest
 infeasible
 maros

Input:
 ($testFileName).txt
Output:
 shown on terminal

".out" files are the standard results to verify the correctness of the program output.
dualtest.txt tests the correctness through comparing the results of simplex with dualsimplex and also compares the running efficiency.