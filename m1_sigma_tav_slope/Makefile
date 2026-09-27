CC=gcc
CFLAGS=-O2
all: m1_slope_sde
m1_slope_sde: m1_slope_sde.c
	$(CC) $(CFLAGS) -o m1_slope_sde m1_slope_sde.c -lm
clean:
	rm -f m1_slope_sde results/*.csv figures/*
