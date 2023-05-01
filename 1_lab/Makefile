build:
	mpic++ cross.cpp -Wall -lm -o cross
run:
	mpirun -np 1 ./cross

build_picture:
	mpic++ cross.cpp -Wall -lm -o cross -DENABLE_SAVING_PICTURE=1
picture: build_picture
	mpirun -np 4 ./cross
	/bin/python3 /home/stanislav/ParallelProgramming/1_lab/script.py

statistic: build
	sudo ./run_lab.sh
	/bin/python3 /home/stanislav/ParallelProgramming/1_lab/statistic.py

clean:
	rm -r cross
clean_picture:
	rm -r res.png res.txt
clean_all: clean clean_picture
