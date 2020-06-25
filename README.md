<p align="left"> 
<img src="https://github.com/d-j-e/SNPPar/blob/master/SNPPar.png" width="150">
</p>

# SNPPar
Parallel/Homoplasic SNP Finder

SNPPar is designed to find homoplasic SNPs based on a user-defined phylogenetic tree - more specifically, it searches for those SNPs that are: parallel - same mutation (eg. A \~> T) @ same position in two (or more) unrelated groups/isolates; convergent - different mutation in resulting in same base (eg. A \~> T, C \~> T) @ same position in two (or more) unrelated groups/isolates; and/or revertant - mutation back to ancestral state (eg. A \~> T \~> A)

By default, SNPPar uses TreeTime for ancestral state reconstruction (ASR), but using FastML for ASR is also available if FastML is installed (though much, much slower)

Current Version: V1.0

# Home:

[SNPPar](https://github.com/d-j-e/SNPPar)

Please use the SNPPar 'Issues' page on GitHub if you have any problems or find any errors in the output. 

Also available: [SNPPar_test](https://github.com/d-j-e/SNPPar_test) - a git with all the data, code, instructions and outputs for testing SNPPar with simulated and empirical data as found in the citation below.

# Citation 
Edwards DJ, Duchêne S, Pope B, and Holt KE. "SNPPar: identifying convergent evolution and other homoplasies from microbial whole-genome alignments" (bioarchive link)

Submitted to _Microbial Genomics_

# License:

[GNU General Public License v3.0](https://github.com/d-j-e/SNPPar/blob/master/LICENSE)

# Language:
[Python3](https://www.python.org/downloads/) (v3.6+)

# Requirements:
* [BioPython](https://biopython.org/) (v1.66+)
* [ETE3](http://etetoolkit.org/)
* [TreeTime](https://github.com/neherlab/treetime)

## Optional Requirement:
[FastML](http://fastml.tau.ac.il/)

# Installing SNPPar

* pip install git+https://github.com/d-j-e/SNPPar

This will also install the requirements above (not FastML).

# Input requirements:

The input tree needs to be bifurcating, rooted (midpoint is fine, but an outgroup is much better...), and in Newick format. Also, the branch lengths in the tree should be substitutions/site not, for example, no. of SNPs on each branch. (Note: you can divide the number of SNPs by the total no. of sites tested to get substitutions/site)

SNPPar takes either a SNP table (a small example is provided below) or a MFASTA file with a second file with the SNP positions (in same order)

SNP table

    Pos,IsolateA,IsolateB,IsolateC
    10,A,A,C
    21,T,C,T
    36,T,T,G
    47,T,-,C

As MFASTA

    >IsolateA
    ATTT
    >IsolateB
    ACT-
    >IsolateC
    CTGC

And a SNP Position file

    10
    21
    36
    47

Note that any ambiguous and missing calls should be indicated by a '-'.

Finally, SNPPar currently requires the GenBank version of the reference genome (same sequence as used to map the reads!). Whilst most annotations include gene tags, if your GenBank reference lacks these, you can either add tags prior to using SNPPar, or (crude) tags will be applied by the program. These will start at "Tag_00001" (number of "0"s depends on number of genes); they are not applied to the GenBank file!

Note: If any gene is split in the reference (including across the origin of the reference), the program will ingnore this gene, but will give a warning to user (see **Logging** below).

# Running SNPPar

    snppar -h
    usage: snppar [-h] [-s SNPTABLE] [-m MFASTA] [-l SNP_POSITION_LIST]
                  [-t TREE] [-g GENBANK] [-E SORTING] [-M MUTATION_EVENTS]
                  [-d DIRECTORY] [-p PREFIX] [-P] [-S] [-C] [-R] [-A] [-a] [-n]
                  [-e] [-u] [-f] [-x FASTML_EXECUTE]
        SNPPar: Parallel/homoplasic SNP Finder V0.4.2dev
    optional arguments:
    -h, --help            show this help message and exit
    -s SNPTABLE, --snptable SNPTABLE
                        SNP table (i.e. RedDog output)
    -m MFASTA, --mfasta MFASTA
                        SNPs in MFASTA format
    -l SNP_POSITION_LIST, --snp_position_list SNP_POSITION_LIST
                        SNP position list (required for MFASTA input)
    -t TREE, --tree TREE  Phylogenetic tree (required)
    -g GENBANK, --genbank GENBANK
                        Genbank reference (required)
    -E SORTING, --sorting SORTING
                        Type of sorting (options: "complex" - slower, less
                        memory; "simple" - faster, more memory; default -
                        intermediate)
    -M MUTATION_EVENTS, --mutation_events MUTATION_EVENTS
                        Mutation events file (previous results)
    -d DIRECTORY, --directory DIRECTORY
                        Output directory
    -p PREFIX, --prefix PREFIX
                        Prefix to add to output files
    -P, --parallel        Flag for reporting of parallel calls
    -S, --strict          Flag to output strict parallel calls (for testing,
                        sets '-P' to True")
    -C, --convergent      Flag for reporting of convergent calls
    -R, --revertant       Flag for reporting of revertant calls
    -A, --all_homoplasic_types
                        Flag for reporting of all three homoplasic types
    -a, --no_all_calls    Flag to turn off reporting of all events at each call
                        position (homoplasic reporting)
    -n, --no_homoplasic   Flag to turn off homoplasic calls output
    -e, --no_all_events   Flag to turn off reporting of all mutation events
    -u, --no_clean_up     Flag to turn off deletion of intermediate files on
                        completion of run
    -f, --fastml          Flag to use fastML for ASR (default ASR: TreeTime)
    -x FASTML_EXECUTE, --fastml_execute FASTML_EXECUTE
                        Command to execute fastML (default command: "fastml"
                        i.e. on PATH)

# Example Commands
## To get homoplasic SNPs including all events reported for each (i.e. default settings):

    snppar -s <alleles.csv> -t <tree.tre> -g <genbank.gbk>

## To only map the SNPs back to the tree:
	
    snppar -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -n 

## To report only parallel events:
	
    snppar -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -P -a -n -e

## To report all three kinds of the homoplasic events:
	
    snppar -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -A

## To get only the homoplasic events (e.g. to remove them)

    snppar -s <alleles.csv> -t <tree.tre> -g <genbank.gbk> -a

# SNPPar sorting
Three versions of the SNP sorting are available when using TreeTime for ASR  

                            Filtered out before ASR  
     
    complex                 singletons and monophyletic SNPs 
                            (tested against tree)  
     
    intermediate (default)  same as complex except non-singleton SNPs
                            with missing calls sent to ASR ()
     
    simple                  singletons only

Complex sorting is the most memory efficient of the three, with simple being about twice as costly (estimate!); intermediate sits somewhere in between (though closer to complex).  
  
Run time is more dependant on missing calls; complex and intermediate sorting are quicker than simple sorting when there are no missing calls. When missing calls are present, complex sorting can be much slower than either simple or intermediate sorting. Intermediate sorting is typically faster than simple.  
  
Complex sorting may be useful when memory is a problem; simple sorting can be used to if you would prefer all the internal SNPs (i.e. non-singletons) to be mapped using ASR.  
  
# SNPPar can also use previous output (mutation events)
A prefix must be added to the run - SNPPar will overwrite results if not careful. Note: SNPPar does not work with the trees it produces, use the original tree

## Example to call the three types of homoplasic SNPs post-run

    snppar -M <all_mutation_events.tsv> -p Run_1_ -t <tree.tre> -P -C -R

# Outputs (default)

* Sequence calls at the internal nodes (MFASTA)
* Mutation tables
  * One with all mutation event calls
  * Another with all mutation events at SNP positions found to be homoplasic
* Tree in NHX (extended Newick) and NEXUS formats
  * Internal node labels (same as found in mutation event tables)
  * Total number of mutation events (SNPs) and homoplasic mutation events on each branch (but see **Important Note** below)
  * NEXUS tree can be read into [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) and [iToL](https://itol.embl.de/)
  * NHX tree can be read by [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html)
* Log file

# Explanation of header in mutation event files
* Common results
  * Position: Position of mutation event in reference seqeunce
  * Type: Intragenic or Intergenic
  * Ancestor_Node: Internal node that is the parent node of the derived node
  * Derived_Node: Node that has mutation - can be internal node or leaf
  * Ancestor_Call: Base call found at the ancestor node
  * Derived_Call: Base call found at derived node - indicates mutation 
* Intragenic
  * Gene: Gene where mutation event is found (identifier: GenBank tag)
  * Strand: Strand which the gene occurs on - 1: Forward Strand, -1: Reverse Strand
  * Codon: Codon in CDS that has mutation
  * Codon_Position: Position within codon that has mutation
  * Ancestor_Codon: Codon found at the ancestor node
  * Derived_Codon: Codon found at the derived node
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
  
# Logging  
SNPPar now includes logging of all (expected) events to a log file. There are three levels of messages; 'INFO' (Information), 'WARNING', and 'CRITICAL'. All three are *always* reported in the log file.  
* 'INFO' are regular runtime messages.  
* 'WARNING' are for problems such as invariant SNP calls or split genes in the GenBank reference which do not affect the running of SNPPar. However, these are excluded in either case, which may affect the user experience(!)  
* 'CRITICAL' are for problems which result in the immediate termination of the program. These need to be resolved before SNPPar will run successfully.  
  
# Test Data  
In the folder 'test_data' is a SNP table and phylogenetic tree from the simulated data set. These, along with the genbank reference 'NC_00962_3_1.gbk', can be used to test your installation. The expected outputs are included in the subfolder 'test_data/test_outputs'.  
  
## Command to run test data  
Navigate from the SNPPar github folder to test_data:  
  
    cd test_data
  
Then to run SNPPar:
  
    snppar –s MTB_Global_L2_alleles.csv -t MTB_Global_L2.tre -g NC_00962_3_1.gbk -d testing

## Example tree from test_data (using [FigTree](http://tree.bio.ed.ac.uk/software/figtree/))
<p align="left"> 
<img src="https://github.com/d-j-e/SNPPar/blob/master/example_node_labelled_nexus.tre.jpg" width="800">
</p>

# Important Note
SNPPar is very accurate, BUT calls where the ancestor is the root node ('N1') are arbituarly assigned (as with any ASR approach). As such, the output trees have no homoplasic events (parallel, convergent, or revertant) mapped to root node, though the total number of SNPs on each branch is estimated using the ratio of the distance to the child nodes of 'N1'.

When a homoplasic event does occur at the root node and is removed, if there is only one other mutation event at the same SNP position, that mutation event is *not* removed from the tree. Keep this in mind when interpreting the tree output.  
