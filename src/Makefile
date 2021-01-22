ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2
  CUDA_FLAGS=-O3
  OMP_FLAGS=-qopenmp
else
  CC=gcc
  CFLAGS=-O3 -mcmodel=medium -march=native -std=gnu99 -Wno-unknown-pragmas
  CUDA_FLAGS=-O3
  OMP_FLAGS=-fopenmp
endif
MPICC=mpicc
NVCC=nvcc
LDFLAGS=-lm

serial: libodp.a
threads: libodp_threads.a
mpi: libodp_mpi.a
mpi_threads: libodp_mpi_threads.a
cuda: libodp_cuda.a
mpi_cuda: libodp_mpi_cuda.a
all: serial threads mpi mpi_threads cuda mpi_cuda

libodp.a: ODP_Set_aspl.o aspl.o utils.o matmul.o
	$(AR) r $@ $^
libodp_threads.a: ODP_Set_aspl.o aspl_threads.o utils.o matmul_threads.o
	$(AR) r $@ $^
libodp_mpi.a: ODP_Set_aspl_mpi.o aspl_mpi.o utils.o matmul.o
	$(AR) r $@ $^
libodp_mpi_threads.a: ODP_Set_aspl_mpi.o aspl_mpi_threads.o utils.o matmul.o
	$(AR) r $@ $^
libodp_cuda.a: ODP_Set_aspl_cuda.o aspl_cuda.o utils_cuda.o utils.o
	$(AR) r $@ $^
libodp_mpi_cuda.a: ODP_Set_aspl_mpi_cuda.o aspl_mpi_cuda.o utils_cuda.o utils.o
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
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
aspl_mpi_cuda.o: aspl_mpi_cuda.cu common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $< -ccbin mpicc
utils.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
matmul.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
matmul_threads.o: utils.c common.h parameter.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ -c $<
utils_cuda.o: utils_cuda.cu common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
ODP_Set_aspl.o: ODP_Set_aspl.c common.h parameter.h
	$(CC) $(CFLAGS) -o $@ -c $<
ODP_Set_aspl_mpi.o: ODP_Set_aspl_mpi.c common.h parameter.h
	$(MPICC) $(CFLAGS) -o $@ -c $<
ODP_Set_aspl_cuda.o: ODP_Set_aspl_cuda.c common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
ODP_Set_aspl_mpi_cuda.o: ODP_Set_aspl_mpi_cuda.c common.h parameter.h
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $< -ccbin mpicc

clean:
	rm -f *.o *.a
