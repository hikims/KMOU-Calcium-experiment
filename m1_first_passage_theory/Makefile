CC=gcc
CFLAGS=-O2
all: m1_fpt_jacobi m1_full_sde
m1_fpt_jacobi: m1_fpt_jacobi.c
	$(CC) $(CFLAGS) -o m1_fpt_jacobi m1_fpt_jacobi.c -lm
m1_full_sde: m1_full_sde.c
	$(CC) $(CFLAGS) -o m1_full_sde m1_full_sde.c -lm
clean:
	rm -f m1_fpt_jacobi m1_full_sde results/*.csv figures/*
