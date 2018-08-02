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

-  R packages: ``cluster``


Other dependencies
~~~~~~~~~~~~~~~~~~~

For using the phylogenetic component, some programs are required and should be in the PATH:
   -   `BLAST <tp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.26/>`_ 2.2.26
   -   `OrthoMcl pipeline <https://github.com/apetkau/orthomcl-pipeline>`_
   -   Package emboss. For Linux (Ubuntu) this is installed through:
::

    apt-get install emboss


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

 
