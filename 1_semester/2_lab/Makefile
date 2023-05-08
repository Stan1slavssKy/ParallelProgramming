SRCDIR = src
ODIR = build

.PHONY: build

build:
	g++ ${SRCDIR}/integration.cpp -Wall -lm -pthread -o ${ODIR}/integration

run:
	${ODIR}/integration 12 1e-12

build_with_result:
	g++ ${SRCDIR}/integration.cpp -Wall -lm -pthread -o ${ODIR}/integration -DPRINT_RESULT_INTEGRAL_ON=1

run_with_result: build_with_result
	${ODIR}/integration 12 1e-12

statistic: build
	sudo ./run_lab.sh
	/bin/python3 /home/stanislav/ParallelProgramming/1_semester/2_lab/statistic.py

clean:
	rm -r ${ODIR}/integration ${ODIR}/delay
