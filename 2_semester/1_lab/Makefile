.PHONY: build

SRCDIR = src
ODIR = build

init_variant:
	g++ ${SRCDIR}/init_variant.cpp -Wall -o ${ODIR}/init_variant

run_init_variant: init_variant
	./${ODIR}/init_variant

baseline_mpi:
	mpic++ ${SRCDIR}/baseline_mpi.cpp -Wall -o ${ODIR}/baseline_mpi

run_baseline_mpi: baseline_mpi
		mpirun -np 10 ${ODIR}/baseline_mpi

statistic_baseline_mpi: baseline_mpi
	./run_lab_mpi.sh ${BASE_TH_NMB} baseline_mpi
	/bin/python3 /home/stanislav/ParallelProgramming/2_semester/1_lab/statistic.py ${BASE_TH_NMB} baseline_mpi

first_task_mpi:
	mpic++ ${SRCDIR}/first_task.cpp -Wall -o ${ODIR}/first_task_mpi

run_first_task_mpi: first_task_mpi
		mpirun -np 8 ${ODIR}/first_task_mpi

statistic_first_task_mpi: first_task_mpi
	./run_lab_mpi.sh ${BASE_TH_NMB} first_task_mpi
	/bin/python3 /home/stanislav/ParallelProgramming/2_semester/1_lab/statistic.py ${BASE_TH_NMB} first_task_mpi

clean:
	rm -r ${ODIR}/*

# clean_statistic:
# 	rm -r ${ODIR}/statistic.txt
# clean_picture_file:
# 	rm -r ${ODIR}/res.txt
# clean_all: clean clean_statistic clean_picture_file
