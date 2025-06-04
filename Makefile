CC = gcc
CFLAGS = -I
DEPS = recon.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

recon: recon.o
	$(CC) -o recon recon.o
