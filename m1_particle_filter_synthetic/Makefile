CC=gcc
CFLAGS=-O2
all: m1_pf_sde
m1_pf_sde: m1_pf_sde.c
	$(CC) $(CFLAGS) -o m1_pf_sde m1_pf_sde.c -lm
clean:
	rm -f m1_pf_sde results/*.csv figures/*
