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

/* scoring scheme */
int score (char a, char b) {
	
	int		score;                              // final score
    int		match, mismatch, gap;               // award scores

	match = 4;
	mismatch = -3;
	gap = -4;

	if ( a == '-' || b == '-' )
		score = gap;
	else
		score = ( a == b ) ? match : mismatch;
	
	return score;

}

/* max3: return the max score */
int max3 (int d, int l, int u) {
	
	int		max;                             // max score

	max = 0;
	max = ( d > max ) ? d : max;
	max = ( l > max ) ? l : max;
	max = ( u > max ) ? u : max;

	return max;
	
} 

/* mss: return the max score source */
char mss (int d, int l, int u) {

	int		max;                             // max score
	char	mss;                             // max score source
	
	max = max3( d, l, u );

	mss = ' ';
	if (max == d) mss = '\\';
	if (max == l) mss = '|';
	if (max == u) mss = '_';

	return mss;

}

void revseq (char seq[]) {

	int		i;                               // index for seq and revseq
	int		l;                               // seq length

	l = strlen( seq );
	
	char	hold[ l ];                       // hold seq
	
	strcpy(hold, seq);

	for ( i = 0; i < l - 1; i++ ) {
		seq[i] = hold[l-1-i]; 
	}

}

/* smith-waterman alignment function */
void smalign (char seq1[], char seq2[], char bts1[], char bts2[]) {

	int		mat[ MAXSEQ + 1 ][ MAXSEQ + 1 ]; // score matrix
	int		i, j;                            // matrix indices
	int		m, n;                            // sequences length
	char	btm[ MAXSEQ + 1 ][ MAXSEQ + 1 ]; // back trace mark
	int		d, u, l;                         // source score
	int		me, mi, mj;                      // max element and its indices
	char	bi;                              // back trace index, forward

	/* initialization the matrix */
	m = strlen( seq1 );
	n = strlen( seq2 );

	for ( i = 0; i <= m; i++ ) mat[i][0] = 0;
	for ( j = 0; j <= n; j++ ) mat[0][j] = 0;

	/* filling the matrix */
	for ( i = 1; i <= m; i ++ )
		for ( j = 1; j <= n; j++ ) {

			// diagonal
			d = mat[i-1][j-1] + score( seq1[i-1], seq2[j-1] ); 
			// up
			u = mat[i-1][j] + score( seq1[i-1], '-' ); 
			// left
			l = mat[i][j-1] + score( '-', seq2[j-1] );
			// find the max score
			mat[i][j] = max3( d, u, l );
			// mark the trace path
			btm[i][j] = mss( d, u, l );

		}

	/* back trace  */

	// print the back trace matrix
	//for ( i = 0; i <= m; i++ ) {
	//	for ( j = 0; j <= n; j++ ) printf("%2d %c\t", mat[i][j], btm[i][j]);
	//	printf("\n");
	//}
	
	// find max element and its indices
	me = 0;
	for ( i = 1; i <= m; i++ )
		for ( j = 1; j <= n; j++  ) {
			if ( mat[i][j] > me ) {
				me = mat[i][j];
				mi = i;
				mj = j;
			}
		}

	// back trace
	for ( i = mi, j = mj, bi = 0 ; i >= 0; bi++ )
		if ( btm[i][j] == '\\' ) {
			bts1[bi] = seq1[i - 1];
			bts2[bi] = seq2[j - 1];
			i--;
			j--;
		} else if ( btm[i][j] == '_')  {
			bts1[bi] = '-';
			bts2[bi] = seq2[j - 1];
			j--;
		} else if ( btm[i][j] == '|' ) {
			bts1[bi] = seq1[i - 1] ;
			bts2[bi] = '-';
			i--;
		} else {
			bts1[bi] = '\0';
			bts2[bi] = '\0';
			revseq(bts1);
			revseq(bts2);
			return;
		}

}

int main ( int argc, char *argv[] ) {

	char	seq1[ MAXSEQ ];                  // input seq1
	char	seq2[ MAXSEQ ];                  // input seq2
	char	bts1[ MAXSEQ ];                  // back trace seq1
	char	bts2[ MAXSEQ ];                  // back trace seq2
   
	if ( argc != 3 ) {

		printf("Error: insufficent arguments\n");
		usage();

	}

	readseq( argv[1], seq1 );
	readseq( argv[2], seq2 );

	smalign( seq1, seq2, bts1, bts2 );

	printf("orignal sequence:\n");
	puts(seq1);
	puts(seq2);
	printf("smith-waterman alignment\n");
	puts(bts1);
	puts(bts2);

}
