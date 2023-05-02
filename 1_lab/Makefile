SRCDIR = src
ODIR = build

.PHONY: build

build:
	mpic++ ${SRCDIR}/cross.cpp -Wall -lm -o ${ODIR}/cross
run:
	mpirun -np 4 ${ODIR}/cross

build_picture:
	mpic++ ${SRCDIR}/cross.cpp -Wall -lm -o ${ODIR}/cross -DENABLE_SAVING_PICTURE=1
picture: build_picture
	mpirun -np 7 ${ODIR}/cross
	/bin/python3 /home/stanislav/ParallelProgramming/1_lab/script.py

statistic: build
	sudo ./run_lab.sh
	/bin/python3 /home/stanislav/ParallelProgramming/1_lab/statistic.py

clean:
	rm -r ${ODIR}/cross
clean_statistic:
	rm -r ${ODIR}/statistic.txt
clean_picture_file:
	rm -r ${ODIR}/res.txt

clean_all: clean clean_statistic clean_picture_file
