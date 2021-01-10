ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -mavx2
  OMP_FLAGS=-qopenmp
  CUDA_BIN=-ccbin icc
else
  CC=gcc
  CFLAGS=-O3
  OMP_FLAGS=-fopenmp
  CUDA_BIN=-ccbin gcc
endif
MPICC=mpicc
NVCC=nvcc
CUDA_FLAGS=-O3
CUDA_MPIBIN=-ccbin mpicc

##############
LIBS=../src/libapsp.a -lm
THREADS_LIBS=../src/libapsp_threads.a -lm
MPI_LIBS=../src/libapsp_mpi.a -lm
MPI_THREADS_LIBS=../src/libapsp_mpi_threads.a -lm
CUDA_LIBS=../src/libapsp_cuda.a -lcuda -lm
MPI_CUDA_LIBS=../src/libapsp_mpi_cuda.a -lcuda -lm
TEST_GENERAL=n12d3.edges n16d4.edges n72d4.edges n100d10.edges n256d5.edges n256d10.edges
TEST_GRID=w4h4d3r2.edges w5h5d4r2.edges w10h10d6r3.edges w10h10d8r4.edges
TMP := $(shell mktemp)
##############
serial:      general.x grid.x non_regular_general.x regular_general.x simple_general.x simple_grid.x random_general.x random_grid.x sa_general.x
threads:     threads_general.x threads_grid.x threads_simple_general.x threads_simple_grid.x
mpi:         mpi_general.x mpi_grid.x  mpi_simple_general.x mpi_simple_grid.x
mpi_threads: mpi_threads_general.x mpi_threads_grid.x mpi_threads_simple_general.x mpi_threads_simple_grid.x
cuda:        cuda_general.x cuda_grid.x cuda_simple_general.x cuda_simple_grid.x
mpi_cuda:    mpi_cuda_general.x mpi_cuda_grid.x mpi_cuda_simple_general.x mpi_cuda_simple_grid.x
all: serial threads mpi mpi_threads cuda mpi_cuda
##############
general.x: general.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
grid.x: grid.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
non_regular_general.x: non_regular_general.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
regular_general.x: regular_general.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
random_general.x: random_general.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
random_grid.x: random_grid.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
simple_general.x: simple_general.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
simple_grid.x: simple_grid.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
sa_general.x: sa_general.c apsp.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
##############
threads_general.x: general.c apsp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_grid.x: grid.c apsp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_simple_general.x: simple_general.c apsp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_simple_grid.x: simple_grid.c apsp.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
##############
mpi_general.x: mpi_general.c apsp.h
	$(MPICC) $(CFLAGS) -o $@ $< $(MPI_LIBS)
mpi_grid.x: mpi_grid.c apsp.h
	$(MPICC) $(CFLAGS) -o $@ $< $(MPI_LIBS)
mpi_simple_general.x: mpi_simple_general.c apsp.h
	$(MPICC) $(CFLAGS) -o $@ $< $(MPI_LIBS)
mpi_simple_grid.x: mpi_simple_grid.c apsp.h
	$(MPICC) $(CFLAGS) -o $@ $< $(MPI_LIBS)
##############
mpi_threads_general.x: mpi_general.c apsp.h
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(MPI_THREADS_LIBS)
mpi_threads_grid.x: mpi_grid.c apsp.h
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(MPI_THREADS_LIBS)
mpi_threads_simple_general.x: mpi_simple_general.c apsp.h
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(MPI_THREADS_LIBS)
mpi_threads_simple_grid.x: mpi_simple_grid.c apsp.h
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(MPI_THREADS_LIBS)
##############
cuda_general.x: cuda_general.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(CUDA_LIBS) $(CUDA_BIN)
cuda_grid.x: cuda_grid.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(CUDA_LIBS) $(CUDA_BIN)
cuda_simple_general.x: cuda_simple_general.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(CUDA_LIBS) $(CUDA_BIN)
cuda_simple_grid.x: cuda_simple_grid.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(CUDA_LIBS) $(CUDA_BIN)
##############
mpi_cuda_general.x: mpi_cuda_general.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(MPI_CUDA_LIBS) $(CUDA_MPIBIN)
mpi_cuda_grid.x: mpi_cuda_grid.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(MPI_CUDA_LIBS) $(CUDA_MPIBIN)
mpi_cuda_simple_general.x: mpi_cuda_simple_general.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(MPI_CUDA_LIBS) $(CUDA_MPIBIN)
mpi_cuda_simple_grid.x: mpi_cuda_simple_grid.c apsp.h
	$(NVCC) $(CUDA_FLAGS) -o $@ $< $(MPI_CUDA_LIBS) $(CUDA_MPIBIN)
##############
test: serial_test
test_all: serial_test threads_test mpi_test mpi_threads_test cuda_test mpi_cuda_test
serial_test: general.x grid.x simple_general.x simple_grid.x
	for i in $(TEST_GENERAL); do ./general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    ./grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	for i in $(TEST_GENERAL); do ./simple_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    ./simple_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	rm -f $(TMP)

threads_test: threads_general.x threads_grid.x threads_simple_general.x threads_simple_grid.x
	for i in $(TEST_GENERAL); do ./threads_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    ./threads_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	for i in $(TEST_GENERAL); do ./threads_simple_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    ./threads_simple_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	rm -f $(TMP)

mpi_test: mpi_general.x mpi_grid.x mpi_simple_general.x mpi_simple_grid.x
	for i in $(TEST_GENERAL); do mpirun -np 2 ./mpi_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    mpirun -np 2 ./mpi_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	for i in $(TEST_GENERAL); do mpirun -np 2 ./mpi_simple_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    mpirun -np 2 ./mpi_simple_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	rm -f $(TMP)

mpi_threads_test: mpi_threads_general.x mpi_threads_grid.x mpi_threads_simple_general.x mpi_threads_simple_grid.x
	for i in $(TEST_GENERAL); do mpirun -np 2 ./mpi_threads_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    mpirun -np 2 ./mpi_threads_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	for i in $(TEST_GENERAL); do mpirun -np 2 ./mpi_threads_simple_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID); do    mpirun -np 2 ./mpi_threads_simple_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	rm -f $(TMP)

cuda_test: cuda_general.x cuda_grid.x cuda_simple_general.x cuda_simple_grid.x
	for i in $(TEST_GENERAL); do ./cuda_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID);    do ./cuda_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	for i in $(TEST_GENERAL); do ./cuda_simple_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID);    do ./cuda_simple_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	rm -f $(TMP)

mpi_cuda_test: mpi_cuda_general.x mpi_cuda_grid.x mpi_cuda_simple_general.x mpi_cuda_simple_grid.x
	for i in $(TEST_GENERAL); do ./mpi_cuda_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID);    do ./mpi_cuda_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	for i in $(TEST_GENERAL); do ./mpi_cuda_simple_general.x ./graphs/general/$$i; done > $(TMP)
	diff $(TMP) ./graphs/general/ASPL.txt
	for i in $(TEST_GRID);    do ./mpi_cuda_simple_grid.x ./graphs/grid/$$i; done > $(TMP)
	diff $(TMP) ./graphs/grid/ASPL.txt
	rm -f $(TMP)

clean:
	rm -f *.x *.o *~
