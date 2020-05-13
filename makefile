CC=gcc
CFLAGS=-Wall -lm

main : main_iplib.c bmp.o ip_lib.o
	$(CC) $^ -o $@ $(CFLAGS)

tests : bmp.o ip_lib.o tests.c
	$(CC) $^ -o $@ $(CFLAGS)

mandelbrot : bmp.o test_bmp.c
	$(CC) $^ -o $@ $(CFLAGS)

lib: bmp.c ip_lib.c
	$(CC) $^ -c $(CFLAGS)

clean: 
	rm bmp.o ip_lib.o mandelbrot tests main_iplib main
