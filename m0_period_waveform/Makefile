CC=gcc
CFLAGS=-O2
all: m0_period_waveform_sweep
m0_period_waveform_sweep: m0_period_waveform_sweep.c
	$(CC) $(CFLAGS) -o m0_period_waveform_sweep m0_period_waveform_sweep.c -lm
clean:
	rm -f m0_period_waveform_sweep results/*.csv figures/*.png
