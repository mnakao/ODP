CC=gcc
CFLAGS=-O3 -I../../include
CLIBS=../../src/libodp.a -lm

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
