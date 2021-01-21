CC=gcc
CFLAGS=-O3
CLIBS=../../src/libodp.a

all: general.x grid.x general_s.x
general.x: general.c common.o odp.h common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
general_s.x: general_s.c common_s.o odp.h common_s.h
	$(CC) $(CFLAGS) common_s.o -o $@ $< $(CLIBS)
grid.x: grid.c common.o odp.h common.h
	$(CC) $(CFLAGS) common.o -o $@ $< $(CLIBS)
common.o: common.c common.h
	$(CC) $(CFLAGS) -c $<
common_s.o: common_s.c common_s.h
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f *.x *.o *~
