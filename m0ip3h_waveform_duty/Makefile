CC=gcc
CFLAGS=-O2
all: m0ip3h_waveform_sweep
m0ip3h_waveform_sweep: m0ip3h_waveform_sweep.c
	$(CC) $(CFLAGS) -o m0ip3h_waveform_sweep m0ip3h_waveform_sweep.c -lm
clean:
	rm -f m0ip3h_waveform_sweep results/*.csv results/*.txt figures/*
