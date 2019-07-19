<p align="left"> 
<img src="https://github.com/d-j-e/SNPPar/blob/master/SNPPar.png" width="150">
</p>

# SNPPar
Parallel SNP Finder

SNPPar is designed to find homoplastic SNPs based on a user-defined phylogenetic tree - more specifically, it searches for those SNPs that occur in parallel (same mutation @ same position in two [or more] unrelated groups/isolates)

By default, SNPPar uses TreeTime for ancestral state reconstruction (ASR), but using FastML for ASR is also available (though much, much slower)

Current Version: V0.0.4

# Home:

https://github.com/d-j-e/SNPPar (currently private)

Please message me (David Edwards) on Slack if you have any problems or find any errors in the output. Or you can use the SNPPar 'Issues' page now this is on GitHub... 

Note that this code is definitely still in development.

Coming Very Soon: SNPPar_test - a git with all the data, code, instructions and outputs for testing SNPPar with simulated data as found in the citation below.

# Citation 
coming soon...

# License:

[GNU General Public License v3.0](https://github.com/d-j-e/SNPPar/blob/master/LICENSE)

# Requirements:
[Python3](https://www.python.org/downloads/), [BioPython](https://biopython.org/), [ETE3](http://etetoolkit.org/), [TreeTime](https://github.com/neherlab/treetime) 

## Optional Requirement:
[FastML](http://fastml.tau.ac.il/)

# Installing BioPython, ETE3 and TreeTime
All are available through 'pip'

* pip install biopython

* pip install ete3

* pip install phylo-treetime

# Input requirements:

The input tree needs to be bifurcating, rooted (midpoint is fine, but an outgroup is better...), and in Newick format.

SNPPar currently only takes SNP tables (a small example is provided below):


    Pos,IsolateA,IsolateB,IsolateC
    10,A,A,C
    20,T,C,T
    36,T,T,G


Also, SNPPar currently requires the GenBank version of the reference genome (same sequence as used to map the reads!)

# Running SNPPar
  
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

## To get parallel SNPs with all SNP reported for each position (i.e. default settings!):

python snppar.py -s <alleles.csv> -t <tree.tre> -g <genbank.gbk>

## To only map the SNPs back to the tree:
	
python snppar.py -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -n 

## To get all of the homoplastic events:
	
python snppar.py -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -R -C -H

## To get a list of only the homoplastic events (e.g. to remove them)

python snppar.py -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -a -n -H 

# Outputs (default)

* Sequence calls at the internal nodes (MFASTA)
* Mutation tables
  * One with all mutation event calls
  * Another with all mutation events at SNP positions found to be parallel
* Tree in NHX (extended Newick) and NEXUS formats
  * Internal node labels (same as found in mutation event tables)
  * Total number of mutation events (SNPs) and parallel mutation events on each branch (but see Important Note below)
  * the NEXUS tree can be read into [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) and [iToL](https://itol.embl.de/)
  * the NHX tree can be read by [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) 

# Explanation of header in mutation event files...
* Common results

  * Position: Position of mutation event in reference seqeunce
  * Type: Intragenic or Intergenic
  * Ancestor_Node: Internal node that is the parent node of the derived node
  * Derived_Node: Node that has mutation - can be internal node or leaf
  * Ancestor_Call: Base found in the ancestor node
  * Derived_Call: Base found in derived node - indicates mutation 
* Intragenic
  * Gene: Gene where mutation event is found (identifier: GenBank tag)
  * Strand: Strand which the gene occurs on - 1: Forward Strand, -1: Reverse Strand
  * Codon: Codon in CDS that has mutation
  * Codon_Position: Position within codon that has mutation
  * Ancestor_Codon: Codon found in the ancestor node
  * Derived_Codon: Codon found in the derived node
  * Ancestor_A.A.: Translated amino acid at the ancestor node
  * Derived_A.A.: Translated amino acid at the derived node
  * Change: With regard to A.A. -> S: synonymous; NS: nonsynonymous; Ambiguous
* Intergenic 
  * Up_Gene: Nearest gene upstream (5') of mutation event
  * Up_Gene_Strand: Strand on which upstream gene occurs (same as Strand)
  * Up_Gene_Distance: Base pair distance from mutation event to upstream gene
  * Down_Gene: Nearest gene downstream (3') of mutation event
  * Down_Gene_Strand: Strand on which downstream gene occurs (same as Strand)
  * Down_Gene_Distance: Base pair distance from mutation event to downstream gene

# Important Note
SNPPar is very accurate (evidence in SNPPar_test very soon!), BUT calls where the ancestor is the root node ('N1') are ***extremely unreliable*** - Indeed the tree has no homoplastic events (parallel, convergent, or revertant) mapped to root node, though the total number of SNPs is estimated using the ratio of the distance to the child nodes of 'N1'.

# Another Important Note
Have just realised that when a homoplastic event with the root node as the ancestor node is removed, if there is only one other mutation event at the same SNP position, that mutation event should be removed from the tree too (At least for parallel events - convergent and revertant mutation events are more complex). Affects both tree formats... **Will be fixed asap!** 
