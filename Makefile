ifeq ($(ENV), intel)
  CC=icc
  MPICC=mpiicc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2
  LDFLAGS=-lm
  OMP_FLAGS=-qopenmp
else ifeq ($(ENV), fugaku)
  CC=fccpx
  MPICC=mpifccpx
  CFLAGS=-Kfast -Nclang
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  MPICC=mpicc
  CFLAGS=-O3 -march=native -std=gnu99 -Wno-unknown-pragmas
  LDFLAGS=-lm
  OMP_FLAGS=-fopenmp
endif
NVCC=nvcc
NVCC_FLAGS=-O3

serial: libodp.a
threads: libodp_threads.a
mpi: libodp_mpi.a
mpi_threads: libodp_mpi_threads.a
cuda: libodp_cuda.a
mpi_cuda: libodp_mpi_cuda.a
all: serial threads mpi mpi_threads cuda mpi_cuda

libodp.a: aspl.o utils.o
	$(AR) r $@ $^
libodp_threads.a: aspl_threads.o utils_threads.o
	$(AR) r $@ $^
libodp_mpi.a: aspl_mpi.o utils.o
	$(AR) r $@ $^
libodp_mpi_threads.a: aspl_mpi_threads.o utils_threads.o
	$(AR) r $@ $^
libodp_cuda.a: aspl_cuda.o utils_cuda.o utils.o
	$(AR) r $@ $^
libodp_mpi_cuda.a: aspl_mpi_cuda.o utils_cuda.o utils.o
	$(AR) r $@ $^

aspl.o: aspl.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
aspl_threads.o: aspl.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
aspl_mpi.o: aspl_mpi.c common.h parameter.h
	$(MPICC) $(CFLAGS) -o $@ -c $<
aspl_mpi_threads.o: aspl_mpi.c common.h parameter.h
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
aspl_cuda.o: aspl_cuda.cu common.h parameter.h
	$(NVCC) $(NVCC_FLAGS) -o $@ -c $<
aspl_mpi_cuda.o: aspl_mpi_cuda.cu common.h parameter.h
	$(NVCC) $(NVCC_FLAGS) -o $@ -c $< -ccbin mpicc
utils.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
utils_threads.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
utils_cuda.o: utils_cuda.cu common.h parameter.h
	$(NVCC) $(NVCC_FLAGS) -o $@ -c $<

clean:
	rm -f *.o *.a
