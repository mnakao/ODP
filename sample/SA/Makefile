ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2 -I../../include
  CLIBS=-lm
  OMP_FLAGS=-qopenmp
  NVCC_BIN=-ccbin icc
else ifeq ($(ENV), fugaku)
  CC=fccpx
  CFLAGS=-Kfast -I../../include -Nclang
  CLIBS=-lm
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  CFLAGS=-O3 -march=native -std=gnu99 -Wno-unknown-pragmas -I../../include
  CLIBS=-lm
  OMP_FLAGS=-fopenmp
  NVCC_BIN=-ccbin gcc
endif
NVCC=nvcc
NVCC_CFLAGS=-O3 -I../../include
############
LIBODP=../../libodp.a
LIBODP_THREADS=../../libodp_threads.a
LIBODP_CUDA=../../libodp_cuda.a -lcuda
all: serial
serial:general.x grid.x general_s.x grid_s.x
threads: threads_general.x threads_grid.x threads_general_s.x threads_grid_s.x
cuda: cuda_general.x cuda_grid.x cuda_general_s.x cuda_grid_s.x
############
general.x: general.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP)
general_s.x: general_s.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP)
grid.x: grid.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP)
grid_s.x: grid_s.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP)
############
threads_general.x: general.c common.o common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_THREADS)
threads_general_s.x: general_s.c common.o common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_THREADS)
threads_grid.x: grid.c common.o common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_THREADS)
threads_grid_s.x: grid_s.c common.o common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_THREADS)
############
cuda_general.x: general.c common.o common.h
	$(NVCC) $(NVCC_CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_CUDA) $(NVCC_BIN)
cuda_general_s.x: general_s.c common.o common.h
	$(NVCC) $(NVCC_CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_CUDA) $(NVCC_BIN)
cuda_grid.x: grid.c common.o common.h
	$(NVCC) $(NVCC_CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_CUDA) $(NVCC_BIN)
cuda_grid_s.x: grid_s.c common.o common.h
	$(NVCC) $(NVCC_CFLAGS) common.o -o $@ $< $(CLIBS) $(LIBODP_CUDA) $(NVCC_BIN)

clean:
	rm -f *.x *.o *~
