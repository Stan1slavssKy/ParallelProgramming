build:
	mpic++ exponent.cpp -Wall -lm -lgmpxx -lgmp -o exponent
run:
	mpirun -np 8 ./exponent 10000 > output.txt