CC=gcc
CFLAGS=-O2
all: m1_isi_sde
m1_isi_sde: m1_isi_sde.c
	$(CC) $(CFLAGS) -o m1_isi_sde m1_isi_sde.c -lm
clean:
	rm -f m1_isi_sde results/*.csv figures/*
