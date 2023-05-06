SRCDIR = src
ODIR = build
DELAY_SRCDIR = delay_time_research

.PHONY: build

build:
	mpic++ ${SRCDIR}/cross.cpp -Wall -lm -o ${ODIR}/cross

run:
	mpirun -np 1 ${ODIR}/cross

build_picture:
	mpic++ ${SRCDIR}/cross.cpp -Wall -lm -o ${ODIR}/cross -DENABLE_SAVING_PICTURE=1

picture: build_picture
	mpirun -np 6 ${ODIR}/cross
	/bin/python3 /home/stanislav/ParallelProgramming/1_semester/1_lab/script.py

statistic: build
	sudo ./run_lab.sh
	/bin/python3 /home/stanislav/ParallelProgramming/1_semester/1_lab/statistic.py

delay_research_run:
	mpic++ ${DELAY_SRCDIR}/delay.cpp -Wall -lm -o ${ODIR}/delay
	mpirun -np 2 ${ODIR}/delay 1000000

clean:
	rm -r ${ODIR}/cross ${ODIR}/delay
clean_statistic:
	rm -r ${ODIR}/statistic.txt
clean_picture_file:
	rm -r ${ODIR}/res.txt
clean_all: clean clean_statistic clean_picture_file
