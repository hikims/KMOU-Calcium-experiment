CC=gcc
CFLAGS=-O2
all: m0ip3h_stimulus_sweep
m0ip3h_stimulus_sweep: m0ip3h_stimulus_sweep.c
	$(CC) $(CFLAGS) -o m0ip3h_stimulus_sweep m0ip3h_stimulus_sweep.c -lm
clean:
	rm -f m0ip3h_stimulus_sweep results/*.csv figures/*
