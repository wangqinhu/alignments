Alignments
==========

Basic nucleotide sequences alignment tools implemented in C for biologist.

Features
---------
- read fasta sequence 
- smith-waterman alignment
- needleman-wunsch alignment
- input/out control

Clone & Compile
---------------

```
git clone https://github.com/wangqinhu/alignments.git
cd alignments
make
```

Usage
-----

- Simith-Waterman Alignment

```
swalign <options>

  -h  help
  -i  input sequence file 1
  -j  input sequence file 2
  -a  input sequence a, directly in command line
  -b  input sequence b, directly in command line
  -s  output the raw alignment
```

- Needleman-Wunsch Alignment

```
nwalign <options>

  -h  help
  -i  input sequence file 1
  -j  input sequence file 2
  -a  input sequence a, directly in command line
  -b  input sequence b, directly in command line
  -s  output the raw alignment
```
