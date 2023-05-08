SRCDIR = src
ODIR = build

.PHONY: build

build:
	g++ ${SRCDIR}/integration.cpp -Wall -lm -pthread -o ${ODIR}/cross

run:
	${ODIR}/cross 1 1e-12

clean:
	rm -r ${ODIR}/cross ${ODIR}/delay

