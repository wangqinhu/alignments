smalign: smalign.c
	gcc -o smalign smalign.c -lm

clean:
	rm -f smalign
