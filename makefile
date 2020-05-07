CC=gcc
CFLAGS=-Wall -lm

main_iplib : main_iplib.c bmp.o ip_lib.o
	$(CC) $^ -o $@ $(CFLAGS)

ip_lib.o : ip_lib.c
	$(CC) $^ -c -o $@ $(CFLAGS)

factorial : bmp.o test_bmp.c
	$(CC) $^ -o $@ $(CFLAGS)

bmp.o: bmp.c
	$(CC) $^ -c -o $@ $(CFLAGS)

clean: 
	rm bmp.o factorial ip_lib.o