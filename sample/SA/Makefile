ifeq ($(ENV), intel)
  CC=icc
  CFLAGS=-O3 -std=gnu99 -Wno-unknown-pragmas -mavx2 -I../../include
  CLIBS=../../libodp.a -lm
else ifeq ($(ENV), fujitsu)
  CC=fccpx
  CFLAGS=-Kfast -D_FUJITSU -I../../include
  CLIBS=../../libodp.a -lm
else
  CC=gcc
  CFLAGS=-O3 -march=native -std=gnu99 -Wno-unknown-pragmas -I../../include
  CLIBS=../../libodp.a -lm
endif

all: general.x grid.x general_s.x grid_s.x general_vt.x
general.x: general.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
general_s.x: general_s.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
general_vt.x: general_vt.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
grid.x: grid.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
grid_s.x: grid_s.c common.o common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
common.o: common.c common.h
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f *.x *.o *~
