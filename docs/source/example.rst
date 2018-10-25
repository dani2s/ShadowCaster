Example
=======

This documentation aims to be a complete example walk through for the usage of ShadowCaster. 
It assumes you have successfully gone through the :doc:`installation`.

Downloading test data
---------------------
Download and extract the test data repository of ShadowCaster in a suitable location, see `here <https://github.com/dani2s/ShadowCaster_testData>`_

Software specifications 
-----------------------

The results provided in the test data repository were obtained running ShadowCaster with the following software versions.

:O.S: Linux Mint 18.03 (Sylvia)
:Package base: Ubuntu Xenial  
:Python: v2.7.15
:R: v3.4.4
:Perl: v5.22.1
:mysql: v5.7.23
:Python packages:
  * ``pip`` v18.0
  * ``numpy`` v1.14.3
  * ``biopython`` v1.72
  * ``matplotlib`` v2.2.2
  * ``pandas`` v0.23.0
  * ``scipy`` v1.1.0
  * ``scikit-learn`` v0.19.1
  * ``ete3`` v3.1.1
:R package:
  * ``cluster`` v2.0.7-1 


Run ShadowCaster
----------------

All fasta files are in the ``shadowcaster-input`` and ``proteomes-output`` folders. The arguments needed in the args.ini file are:

  * query_genome = /home/user/path/to/shadowcaster-input/Rdenitrificans_genome.fasta
  * query_proteome = /home/user/path/to/shadowcaster-input/Rhodanobacter_denitrificans.fasta
  * proteomes_folder = /home/user/path/to/proteomes-output/proteomes/
  * orthomcl_config = /home/user/orthomcl-pipeline/orthomcl.conf
  * blastp26 = Binary file or your command line used of blastp v2.2.26.
  * formatdb26 = Binary file or your command line used of blastall v2.2.26.
  * nuSVM = 0.4

MUST use the full path of the files or directory.

Run ShadowCaster through: 

::

    cd $ShadowCaster/bin/ 
    shadowcaster --config_file args.ini


When ShadowCaster has finished the message "ShadowCaster finished" 
is printed. The program generates a number of
files in the output directory (called with the date and time of the running).

Output description
------------------

ShadowCaster generates the following files in the output directory.

:log.txt: Contains used parameters

:Parametric folder:
  * ``data4mer_chi2.csv`` Compositional difference measured with chi2 between each gene and the genome based on 4mers . 
  * ``kl_codonUsage.csv`` Measure compositional difference of codon usage with Kullback-Leibler.
  * ``plot_aliens.png`` Plot of the two compositional features.
  * ``alien_genes.fasta`` Fasta file with the alien genes identified by one-class support vector machine.
  
:Phylogenetic folder:
  * Orthomcl folder contains files of orthologs groups found by OrthoMCL.
  * ``orthologs_probabilities.csv`` the probabilities of orthology between the query proteome and each other genome in the phylogenetic shadow.
  * ``alien_likelihoods.csv`` The log likelihood calculated for each alien gene found in the parametric component.
  * ``histogram_alienLikelihoods.png`` Histogram of likelihoods from alien genes.
  * ``hgt_candidates.csv`` List of HGT candidate genes found with the software with their likelihood.
  * ``shadowcaster_predictions.fasta`` Fasta file with genes predicted as HGT by ShadowCaster.


A copy of all the output files can be found in the ``shadowcaster-output`` folder.
