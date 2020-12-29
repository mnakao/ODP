ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas
  CUDA_FLAGS=-O3
  OMP_FLAGS=-qopenmp
else
  CC=gcc
  CFLAGS=-O3 -mcmodel=medium -march=native -std=gnu99
  CUDA_FLAGS=-O3
  OMP_FLAGS=-fopenmp
endif
MPICC=mpicc
NVCC=nvcc
LDFLAGS=-lm

serial: libapsp.a
threads: libapsp_threads.a
mpi: libapsp_mpi.a
mpi_threads: libapsp_mpi_threads.a
cuda: libapsp_cuda.a
mpi_cuda: libapsp_mpi_cuda.a
all: serial threads mpi mpi_threads cuda mpi_cuda

libapsp.a: all_run.o apsp.o utils.o
	$(AR) r $@ $^
libapsp_threads.a: all_run.o apsp_threads.o utils_threads.o
	$(AR) r $@ $^
libapsp_mpi.a: all_mpi_run.o apsp_mpi.o utils.o
	$(AR) r $@ $^
libapsp_mpi_threads.a: all_mpi_run.o apsp_mpi_threads.o utils_threads.o
	$(AR) r $@ $^
libapsp_cuda.a: all_cuda_run.o apsp_cuda.o utils_cuda.o utils.o
	$(AR) r $@ $^
libapsp_mpi_cuda.a: all_mpi_cuda_run.o apsp_mpi_cuda.o utils_cuda.o utils.o
	$(AR) r $@ $^

apsp.o: apsp.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
apsp_threads.o: apsp.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
apsp_mpi.o: apsp_mpi.c common.h parameter.h
	$(MPICC) $(CFLAGS) -o $@ -c $<
apsp_mpi_threads.o: apsp_mpi.c common.h parameter.h
	$(MPICC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
apsp_cuda.o: apsp_cuda.cu common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
apsp_mpi_cuda.o: apsp_mpi_cuda.cu common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $< -ccbin mpicc
utils.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
utils_threads.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
utils_cuda.o: utils_cuda.cu common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
all_run.o: all_run.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
all_mpi_run.o: all_mpi_run.c common.h parameter.h
	$(MPICC) $(CFLAGS) -o $@ -c $<
all_cuda_run.o: all_cuda_run.c common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
all_mpi_cuda_run.o: all_mpi_cuda_run.c common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $< -ccbin mpicc

clean:
	rm -f *.o *.a
