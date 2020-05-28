CC=gcc
CFLAGS=-Wall -lm
CFLAGS_MAIN=-Wall -ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra



main : src/main_iplib.c lib/bmp.o lib/ip_lib.o 
	$(CC) $^ -o $@ $(CFLAGS_MAIN)

lib/ip_lib.o : src/ip_lib.c
	$(CC) $^ -c -o $@ $(CFLAGS)

lib/bmp.o : src/bmp.c
	$(CC) $^ -c -o $@ $(CFLAGS)

test_mat : target/bmp.o target/ip_lib.o test/test_mat.c
	$(CC) $^ -o test/$@ $(CFLAGS)

mandelbrot : lib/bmp.o lib/ip_lib.o test/test_bmp.c
	$(CC) $^ -o test/$@ $(CFLAGS)


clean: 
	rm lib/* main
	rm test/test_bmp test/test_mat
