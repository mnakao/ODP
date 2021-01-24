CC=gcc
CFLAGS=-O3
CLIBS=../../src/libodp.a

all: general.x grid.x general_s.x
general.x: general.c common.o odp.h common.h utils.o
	$(CC) $(CFLAGS) common.o utils.o -o $@ $< $(CLIBS)
general_s.x: general_s.c common.o odp.h common.h utils.o
	$(CC) $(CFLAGS) common.o utils.o -o $@ $< $(CLIBS)
grid.x: grid.c common.o odp.h common.h utils.o
	$(CC) $(CFLAGS) common.o utils.o -o $@ $< $(CLIBS)
common.o: common.c common.h
	$(CC) $(CFLAGS) -c $<
utils.o: utils.c
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f *.x *.o *~
