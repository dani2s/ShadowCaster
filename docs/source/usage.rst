Usage
=====

ShadowCaster uses a configuration file (args.ini) to manage all the options needed.
These file must be specified on the command line and will supply the following arguments:

Files
-----
- **Query genome** *Only fasta files*
- **Query proteome** *Only fasta files*
- **OrthoMcl configuration file**

Path
----
- **Proteomes** folder to calculate the phylogenetic shadow *Only fasta files*
- **Blastp26** *ShadowCaster uses blastp 2.2.26, specify the binary file or the command line used*
- **Formatdb** *specify the binary file or the command used in the shell. Ex. formatdb*

Parametric
----------
- **nuSVM** *A bound between the fraction of training errors and the fraction of support vectors. Should be in the interval (0, 1].* 


A template of the args.ini file can be find in the bin folder.