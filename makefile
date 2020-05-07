CC=gcc
CFLAGS=-Wall -lm

factorial : bmp.o test_bmp.c
	$(CC) $^ -o $@ $(CFLAGS)

bmp.o: bmp.c
	$(CC) -c $^ -o $@

clean: 
	rm bmp.o factorial