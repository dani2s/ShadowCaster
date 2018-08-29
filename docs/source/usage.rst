Usage
=====

ShadowCaster uses a configuration file (args.ini) to manage all the options needed.
This file must be specified on the command line and will supply the following arguments:

Files
-----
- **Query genome** *Only fasta files.*
- **Query proteome** *Only fasta files*
- **OrthoMcl configuration file**, previously obtained with the installation of OrthoMcl-pipeline(orthomcl.conf).

Path
----
- **Proteomes folder**, contain the proteomes(*Only fasta files*) of each species to construct the shadow. 
	-Provided by the user (a collection of FASTA files)
	
	or
	
	-Automatically retrieved by ShadowCaster from the NCBI ftp site by using ``script/get_proteomes.py``, see :doc:`scripts`.


- **Blastp26** ShadowCaster uses blastp 2.2.26, specify the binary file or the shell command used
- **Formatdb** Specify the binary file or the shell command used. Ex. formatdb

Parametric
----------
- **nuSVM** A bound between the fraction of training errors and the fraction of support vectors. Should be in the interval (0, 1]. 


A template of the args.ini file can be found in the bin folder.

Specifications of fasta files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The id of each gene in the GENOME fasta file SHOULD not contain any character like ``|, #, %, /, \, *, &, $, !, :``.
It is preferred only one identification number.

All the proteomes files provided by the user MUST follow these specifications:
  * Fasta file name ONLY can be the binomial name of the species (Rhodanobacter_denitrificans.fasta).
  * Each id record of a fasta file MUST have only one id number or the following structure:

::

    >AGG91012.1 ribosomal protein L34 [Rhodanobacter denitrificans]
    MKRTFQPSKLKRARTHGFRARMATADGRKVLNARRAKGRKRLIP 	


Examples of the input files can be found in the test data repository of ShadowCaster, see `here <https://github.com/dani2s/ShadowCaster_testData>`_

Run ShadowCaster
----------------

::

	cd $ShadowCaster/bin/
	shadowcaster --config_file args.ini
