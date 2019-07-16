### SNPPar
Parallel SNP Finder

V0.0.2

## Home:

https://github.com/d-j-e/SNPPar (currently private)

## Requirements:

Python 3, BioPython, ete3, TreeTime 

# Optional Requirement:

FastML

## Installing biopython, ete3 and TreeTime

All are available through 'pip'

pip install biopython
pip install ete3
pip install phylo-treetime


# To get the options for SNPPar (or see below):

python snppar.py -h


# To get parallel SNPs with all SNP reported for each position (i.e. default settings!):

python snppar.py -s <alleles.csv> -t <tree> -g <genbank.gbk>


# To only map the SNPs back to the tree:
	
python snppar.py -s <alleles.csv> -t <tree> -g <genbank.gbk> -n 


# To get all of the homoplastic events (and any other change(s) at the same positions):
	
python snppar.py -s <alleles.csv> -t <tree> -g <genbank.gbk> -R -C -H


# To get a list of only the homoplastic events (e.g. to remove them)

python snppar.py -s <alleles.csv> -t <tree> -g <genbank.gbk> -a -n -H 


* Please message me (David Edwards) on Slack if you have any problems or find any errors in the output. This code is definitely still in development.

# Warning, the name may yet change...

## Input requirements:

The tree needs to be bifucating, rooted (midpoint is fine, but an outgroup is better...), and in Newick format.

SNPPar currently only takes SNP tables (a small example is provided below):


Pos,IsolateA,IsolateB,IsolateC
10,A,A,C
20,T,C,T
36,T,T,G


Also, SNPPar currently requires the GenBank version of the reference genome (same sequence as used to map the reads!)


python snppar.py -h

usage: snppar.py [-h] -s SNPTABLE -t TREE -g GENBANK [-d DIRECTORY]
                 [-p PREFIX] [-S] [-H] [-C] [-R] [-a] [-n] [-e] [-f] [-c] [-u]

SNPPar: Parallel SNP Finder V0.0.2

optional arguments:
  -h, --help            show this help message and exit
  -s SNPTABLE, --snptable SNPTABLE
                        SNP table (required)
  -t TREE, --tree TREE  Phylogenetic tree (required)
  -g GENBANK, --genbank GENBANK
                        Genbank reference (required)
  -d DIRECTORY, --directory DIRECTORY
                        Output directory
  -p PREFIX, --prefix PREFIX
                        Prefix to add to output files
  -S, --strict          Flag to output strict parallel calls (for testing)
  -H, --homoplastic     Flag for reporting of all homoplastic calls
  -C, --convergent      Flag for reporting of convergent calls
  -R, --revertant       Flag for reporting of revertant calls
  -a, --no_all_calls    Flag to turn off reporting of all mutation events at
                        each call position
  -n, --no_parallel     Flag to turn off parallel calls output
  -e, --no_all_events   Flag to turn off reporting of all mutation events
  -f, --fastml          Flag to use fastML for ASR (default ASR: TreeTime)
  -c, --counting        Flag to display counts during SNP testing
  -u, --no_clean_up     Flag to turn off deletion of intermediate files on
                        completion of run
