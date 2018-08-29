======================
get_proteomes.py
======================

Prerequisites
-------------

* `EDirect <https://www.ncbi.nlm.nih.gov/books/NBK179288/>`_ UNIX command line of NCBI.

Before using the script, check that the commands esearch and xtract work correctly in a new shell window.

::
			
	type esearch xtract
	

Usage
-----
The usage and help documentation of ``get_proteomes.py`` can be seen by
running ``python get_proteomes.py -h``:


Example
-------
An example of how to run ``get_proteomes.py`` on the test data::
    
    cd ShadowCaster/scripts
    python get_proteomes.py -n Rhodanobacter_denitrificans -sp 30

This results in the following output files in the folder named with the species name provided:

    * ``log.txt`` Name of the downloaded species and its ftp address.
    * ``proteomes folder`` Proteomes (fasta file) used to construct the shadow.


The results should be similar to those found in the ``ftp-output`` folder of the test data repository, see `here <https://github.com/dani2s/ShadowCaster_testData>`_

