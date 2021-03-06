ifeq ($(ENV), intel)
  CC=icc
  MPICC=mpiicc
  CFLAGS=-O3 -std=gnu99 -mavx2 -I../include
  OMP_FLAGS=-qopenmp
  NVCC_BIN=-ccbin icc
  NVCC_MPIBIN=-ccbin mpiicc
else ifeq ($(ENV), fugaku)
  CC=fccpx
  MPICC=mpifccpx
  CFLAGS=-Kfast -I../include -Nclang
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  MPICC=mpicc
  CFLAGS=-O3 -I../include
  OMP_FLAGS=-fopenmp
  NVCC_BIN=-ccbin gcc
  NVCC_MPIBIN=-ccbin mpicc
endif
NVCC=nvcc
NVCC_CFLAGS=-O3 -I../include

##############
LIBS=../libodp.a -lm
THREADS_LIBS=../libodp_threads.a -lm
MPI_LIBS=../libodp_mpi.a -lm
MPI_THREADS_LIBS=../libodp_mpi_threads.a -lm
CUDA_LIBS=../libodp_cuda.a -lcuda -lm
MPI_CUDA_LIBS=../libodp_mpi_cuda.a -lcuda -lm
TEST_GENERAL=n12d3.edges n16d4.edges n72d4.edges n100d10.edges n256d5.edges n256d10.edges
TEST_GRID=w4h4d3r2.edges w5h5d4r2.edges w10h10d6r3.edges w10h10d8r4.edges
TMP := $(shell mktemp)
##############
serial:      general.x grid.x non_regular_general.x regular_general.x random_general.x random_grid.x
threads:     threads_general.x threads_grid.x
mpi:         mpi_general.x mpi_grid.x
mpi_threads: mpi_threads_general.x mpi_threads_grid.x
cuda:        cuda_general.x cuda_grid.x
mpi_cuda:    mpi_cuda_general.x mpi_cuda_grid.x
all: serial threads mpi mpi_threads cuda mpi_cuda
##############
general.x: general.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
grid.x: grid.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
non_regular_general.x: non_regular_general.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
regular_general.x: regular_general.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
random_general.x: random_general.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
random_grid.x: random_grid.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
##############
threads_general.x: general.c
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
threads_grid.x: grid.c
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(THREADS_LIBS)
##############
mpi_general.x: mpi_general.c
	$(MPICC) $(CFLAGS) -o $@ $< $(MPI_LIBS)
mpi_grid.x: mpi_grid.c
	$(MPICC) $(CFLAGS) -o $@ $< $(MPI_LIBS)
##############
mpi_threads_general.x: mpi_general.c
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(MPI_THREADS_LIBS)
mpi_threads_grid.x: mpi_grid.c
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(MPI_THREADS_LIBS)
##############
cuda_general.x: general.c
	$(NVCC) $(NVCC_CFLAGS) -o $@ $< $(CUDA_LIBS) $(NVCC_BIN)
cuda_grid.x: grid.c
	$(NVCC) $(NVCC_CFLAGS) -o $@ $< $(CUDA_LIBS) $(NVCC_BIN)
##############
mpi_cuda_general.x: mpi_general.c
	$(NVCC) $(NVCC_CFLAGS) -o $@ $< $(MPI_CUDA_LIBS) $(NVCC_MPIBIN)
mpi_cuda_grid.x: mpi_grid.c
	$(NVCC) $(NVCC_CFLAGS) -o $@ $< $(MPI_CUDA_LIBS) $(NVCC_MPIBIN)

##############
test: serial_test
test_all: serial_test threads_test mpi_test mpi_threads_test cuda_test mpi_cuda_test
serial_test: general.x grid.x
	for i in $(TEST_GENERAL); do ./general.x ./graph/general/$$i; done > $(TMP)
	diff $(TMP) ./graph/general/ASPL.txt
	for i in $(TEST_GRID); do ./grid.x ./graph/grid/$$i; done > $(TMP)
	diff $(TMP) ./graph/grid/ASPL.txt
	rm -f $(TMP)

threads_test: threads_general.x threads_grid.x
	for i in $(TEST_GENERAL); do ./threads_general.x ./graph/general/$$i; done > $(TMP)
	diff $(TMP) ./graph/general/ASPL.txt
	for i in $(TEST_GRID); do ./threads_grid.x ./graph/grid/$$i; done > $(TMP)
	diff $(TMP) ./graph/grid/ASPL.txt
	rm -f $(TMP)

mpi_test: mpi_general.x mpi_grid.x
	for i in $(TEST_GENERAL); do mpirun -np 2 ./mpi_general.x ./graph/general/$$i; done > $(TMP)
	diff $(TMP) ./graph/general/ASPL.txt
	for i in $(TEST_GRID); do mpirun -np 2 ./mpi_grid.x ./graph/grid/$$i; done > $(TMP)
	diff $(TMP) ./graph/grid/ASPL.txt
	rm -f $(TMP)

mpi_threads_test: mpi_threads_general.x mpi_threads_grid.x
	for i in $(TEST_GENERAL); do mpirun -np 2 ./mpi_threads_general.x ./graph/general/$$i; done > $(TMP)
	diff $(TMP) ./graph/general/ASPL.txt
	for i in $(TEST_GRID); do mpirun -np 2 ./mpi_threads_grid.x ./graph/grid/$$i; done > $(TMP)
	diff $(TMP) ./graph/grid/ASPL.txt
	rm -f $(TMP)

cuda_test: cuda_general.x cuda_grid.x
	for i in $(TEST_GENERAL); do ./cuda_general.x ./graph/general/$$i; done > $(TMP)
	diff $(TMP) ./graph/general/ASPL.txt
	for i in $(TEST_GRID); do ./cuda_grid.x ./graph/grid/$$i; done > $(TMP)
	diff $(TMP) ./graph/grid/ASPL.txt
	rm -f $(TMP)

mpi_cuda_test: mpi_cuda_general.x mpi_cuda_grid.x
	for i in $(TEST_GENERAL); do ./mpi_cuda_general.x ./graph/general/$$i; done > $(TMP)
	diff $(TMP) ./graph/general/ASPL.txt
	for i in $(TEST_GRID); do ./mpi_cuda_grid.x ./graph/grid/$$i; done > $(TMP)
	diff $(TMP) ./graph/grid/ASPL.txt
	rm -f $(TMP)

clean:
	rm -f *.x *.o *~
