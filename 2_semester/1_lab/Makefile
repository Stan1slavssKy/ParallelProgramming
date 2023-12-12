.PHONY: build

SRCDIR = src
ODIR = build
PY_SCRIPT_PATH = /home/stanislav/ParallelProgramming/2_semester/1_lab/

init_variant:
	g++ ${SRCDIR}/init_variant.cpp -Wall -o ${ODIR}/init_variant

run_init_variant: init_variant
	./${ODIR}/init_variant

# ==================================================================================================================

baseline_mpi:
	mpic++ ${SRCDIR}/baseline_mpi.cpp -Wall -o ${ODIR}/baseline_mpi

run_baseline_mpi: baseline_mpi
		mpirun -np 10 ${ODIR}/baseline_mpi

statistic_baseline_mpi: baseline_mpi
	./run_lab_mpi.sh ${BASE_TH_NMB} baseline_mpi
	/bin/python3 ${PY_SCRIPT_PATH}/statistic.py ${BASE_TH_NMB} baseline_mpi

# ==================================================================================================================

baseline_openmp:
	g++ ${SRCDIR}/baseline_openmp.cpp -Wall -fopenmp -o ${ODIR}/baseline_openmp

run_baseline_openmp: baseline_openmp
		./${ODIR}/baseline_openmp 8

statistic_baseline_openmp: baseline_openmp
	./run_lab_openmp.sh ${BASE_TH_NMB} baseline_openmp
	/bin/python3 ${PY_SCRIPT_PATH}/statistic.py ${BASE_TH_NMB} statistic_baseline_openmp.txt baseline_openmp

# ==================================================================================================================

first_task_mpi:
	mpic++ ${SRCDIR}/first_task.cpp -Wall -o ${ODIR}/first_task_mpi

run_first_task_mpi: first_task_mpi
		mpirun -np 1 ${ODIR}/first_task_mpi

statistic_first_task_mpi: first_task_mpi
	./run_lab_mpi.sh ${BASE_TH_NMB} first_task_mpi
	/bin/python3 ${PY_SCRIPT_PATH}/statistic.py ${BASE_TH_NMB} statistic_first_task_mpi first_task_mpi

# ==================================================================================================================

second_task_openmp:
	g++ ${SRCDIR}/second_task.cpp -Wall -fopenmp -o ${ODIR}/second_task_openmp

run_second_task_openmp: second_task_openmp
		./${ODIR}/second_task_openmp 10

statistic_second_task_openmp: second_task_openmp
	./run_lab_openmp.sh ${BASE_TH_NMB} second_task_openmp
	/bin/python3 ${PY_SCRIPT_PATH}/statistic.py ${BASE_TH_NMB} statistic_second_task_openmp.txt second_task_openmp

# ==================================================================================================================

third_task_mpi:
	mpic++ ${SRCDIR}/third_task.cpp -Wall -o ${ODIR}/third_task_mpi

run_third_task_mpi: third_task_mpi
		mpirun -np 8 ${ODIR}/third_task_mpi

statistic_third_task_mpi: third_task_mpi
	./run_lab_mpi.sh ${BASE_TH_NMB} third_task_mpi
	/bin/python3 ${PY_SCRIPT_PATH}/statistic.py ${BASE_TH_NMB} statistic_third_task_mpi third_task_mpi

# ==================================================================================================================

clean:
	rm -r ${ODIR}/*

# clean_statistic:
# 	rm -r ${ODIR}/statistic.txt
# clean_picture_file:
# 	rm -r ${ODIR}/res.txt
# clean_all: clean clean_statistic clean_picture_file
