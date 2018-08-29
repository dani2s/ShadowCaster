Installation
============

Prerequisites 
--------------

Fundamental prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~

::

    python v2.7.*
    R => v3.4
    perl => v5.10

These items are prerequisites for the installation of ShadowCaster as
described below. 

Packages
~~~~~~~~

-  Install python packages using:

::

    pip install numpy scipy biopython pandas ete3 scikit-learn matplotlib

-  R packages: ``cluster`` v2.0.*
 
-  Install figlet. For Linux (Mint 18.3) use:

::

    apt install figlet


Other dependencies
~~~~~~~~~~~~~~~~~~~

For using the phylogenetic component, some programs are required and should be in the PATH:

-   `BLAST <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.26/>`_ 2.2.26+ (blastp needed)
   
-   `OrthoMcl pipeline <https://github.com/apetkau/orthomcl-pipeline>`_ 
   
OrthoMcl needs blast v2.2.26(blastall and formatdb), these binary files can be found `here <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/>`_. Download blast-2.2.26-*.tar.gz.

For Linux (Mint 18.3) this can be installed through:

::
			
	apt-get install blast2
     		
   		
-   Package emboss. 

For Linux (Mint 18.3) this is installed through:

::
			
	apt-get install emboss


To check that these dependencies were correctly added in the path, run:

::
			
	type orthomcl-pipeline blastall formatdb
	

Installation of ShadowCaster
----------------------------

To use ShadowCaster, download it from the GitHub repository and extract the
files. If you have git installed, you can install ShadowCaster by running:
::

    cd
    git clone https://github.com/dani2s/ShadowCaster.git
 
Resolve all dependencies, see above and then execute:
::

    cd ShadowCaster/ 
    python setup.py install

This will install ShadowCaster under your home folder.
 
