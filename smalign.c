/*  Smith-Waterman alignment
 *
 *  Wang, Qinhu
 *  Northwest A&F University
 *  wangqinhu@nwafu.edu.cn
 *
 *  Feb 28, 2014
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCOL	100
#define MAXSEQ  1000

void usage() {
	
	printf("usage:\n\n");
	printf("./smalign fasta1 fasta2\n\n");
	exit(1);

}

/* define nucleotude */
int isnuc (char nuc) {

	switch (nuc) {
		case 'A': return 1;
		case 'C': return 1;
		case 'G': return 1;
		case 'T': return 1;
		case 'U': return 1;
		case 'a': return 1;
		case 'c': return 1;
		case 'g': return 1;
		case 't': return 1;
		case 'u': return 1;
		default:  return 0;
	}

}

/* read a fasta file and pass the sequence to an array */
void readseq (char filename[], char seq[]) {

	FILE	*file;                           // hold fasta file
	int		i;                               // index for seq[]
	char	line[ MAXCOL ];                  // store each line
	int		j;                               // index for line[]

	file = fopen( filename, "r");

	if ( file == NULL ) {

		printf("Error: file %s does not exist!\n", filename);
		exit(1);

	}

	i = 0;
	while ( fgets( line, sizeof line, file ) ) {

		/* if a fasta header line was found */
		if ( line[0] == '>' ) continue;
		
		/* now try to parse the fasta body */
		for ( j = 0; j < strlen(line); j++ ) {
			if ( isnuc( line[j] ) ) seq[i++] = line[j];
		}

	} // end of while
	seq[i] = '\0';

	fclose( file );

} // end readseq

/* smith-waterman alignment function */
void smalign (char seq1[], char seq2[]) {

}

int main ( int argc, char *argv[] ) {

	char	seq1[ MAXSEQ ];
	char	seq2[ MAXSEQ ];
   
	if ( argc != 3 ) {

		printf("Error: insufficent argments\n");
		usage();

	}

	readseq( argv[1], seq1 );
	readseq( argv[2], seq2 );

	smalign( seq1, seq2 );

}
