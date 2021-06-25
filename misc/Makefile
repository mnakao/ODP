ifeq ($(ENV), intel)
  CC=icc
  MPICC=mpiicc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2 -I../include
  CLIBS=-lm
  OMP_FLAGS=-qopenmp
  NVCC_BIN=-ccbin icc
  NVCC_MPIBIN=-ccbin mpiicc
else ifeq ($(ENV), fugaku)
  CC=fccpx
  MPICC=mpifccpx
  CFLAGS=-Kfast -I../include -Nclang
  CLIBS=-lm
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  MPICC=mpicc
  CFLAGS=-O3 -march=native -std=gnu99 -Wno-unknown-pragmas -I../include
  CLIBS=-lm
  OMP_FLAGS=-fopenmp
  NVCC_BIN=-ccbin gcc
  NVCC_MPIBIN=-ccbin mpicxx
endif
NVCC=nvcc
NVCC_CFLAGS=-O3 -I../include
############
LIBODP=../libodp.a
LIBODP_MPI=../libodp_mpi.a -lm
LIBODP_THREADS=../libodp_threads.a
LIBODP_CUDA=../libodp_cuda.a -lcuda
LIBODP_MPI_CUDA=../libodp_mpi_cuda.a -lcuda -lm
all: serial
serial:general_3opt.x general_s_2opt.x
threads: threads_general_3opt.x
mpi: mpi_general_3opt.x
cuda: cuda_general_3opt.x
mpi_cuda: mpi_cuda_general_3opt.x
############
general_3opt.x: general_3opt.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP)
general_s_2opt.x: general_s_2opt.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP)
############
threads_general_3opt.x: general_3opt.c common.o common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_THREADS)
##############
mpi_general_3opt.x: mpi_general_3opt.c common.o common.h
	$(MPICC) $(CFLAGS) common.o -o $@ $< $(LIBODP_MPI)
############
cuda_general_3opt.x: general_3opt.c common.o common.h
	$(NVCC) $(NVCC_CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_CUDA) $(NVCC_BIN)
############
clean:
	rm -f *.x *.o *~
