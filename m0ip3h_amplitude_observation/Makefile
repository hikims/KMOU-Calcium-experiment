CC=gcc
CFLAGS=-O2
all: m0ip3h_amplitude_sweep
m0ip3h_amplitude_sweep: m0ip3h_amplitude_sweep.c
	$(CC) $(CFLAGS) -o m0ip3h_amplitude_sweep m0ip3h_amplitude_sweep.c -lm
clean:
	rm -f m0ip3h_amplitude_sweep results/*.csv figures/*
