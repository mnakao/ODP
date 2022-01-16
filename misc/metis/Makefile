ifeq ($(ENV), intel)
  CC=icc
  MPICC=mpiicc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2 -I../../include
  CLIBS=-lm
  OMP_FLAGS=-qopenmp
  NVCC_BIN=-ccbin icc
  NVCC_MPIBIN=-ccbin mpiicc
else ifeq ($(ENV), fugaku)
  CC=fccpx
  MPICC=mpifccpx
  CFLAGS=-Kfast -I../../include -Nclang
  CLIBS=-lm
  OMP_FLAGS=-Kopenmp
else
  CC=gcc
  MPICC=mpicc
  CFLAGS=-O3 -march=native -std=gnu99 -Wno-unknown-pragmas -I../../include
  CLIBS=-lm
  OMP_FLAGS=-fopenmp
  NVCC_BIN=-ccbin gcc
  NVCC_MPIBIN=-ccbin mpicxx
endif
NVCC=nvcc
NVCC_CFLAGS=-O3 -I../../include
############
LIBODP=../../libodp.a
LIBODP_THREADS=../../libodp_threads.a
METIS_LIB=/home/mnakao/work/metis-5.1.0/lib/
METIS_INC=/home/mnakao/work/metis-5.1.0/include/
all: serial
serial:general.x general_s.x
############
general.x: general.c detect_temp.o common.o common.h
	$(CC) $(CFLAGS) detect_temp.o common.o -o $@ $< $(CLIBS) $(LIBODP) -I$(METIS_INC) -L$(METIS_LIB) -lmetis
general_s.x: general_s.c detect_temp.o common.o common.h
	$(CC) $(CFLAGS) detect_temp.o common.o -o $@ $< $(CLIBS) $(LIBODP) -I$(METIS_INC) -L$(METIS_LIB) -lmetis
detect_temp.o: detect_temp.c common.h
	$(CC) $(CFLAGS) -c $< $(CLIBS) -I$(METIS_INC)
threads_general.x: general.c detect_temp.o common.o common.h
	$(CC) $(CFLAGS) $(OMP_FLAGS) detect_temp.o common.o -o $@ $< $(CLIBS) $(LIBODP_THREADS)
############
clean:
	rm -f *.x *.o *~
