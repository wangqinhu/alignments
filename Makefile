all: swalign.c nwalign.c
	gcc -o swalign swalign.c -lm
	gcc -o nwalign nwalign.c -lm

clean:
	rm -f swalign nwalign
