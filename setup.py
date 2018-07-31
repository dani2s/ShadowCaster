#!/usr/bin/env python

from setuptools import setup

version = '0.9-beta'

description = """ShadowCaster is a software that combines two types of 
      information - sequence composition and phylogenetic likelihoods. """

setup(name='ShadowCaster',
      version= version,
      description= "Detects Horizontal Gene Transfer events in prokaryotes",
      long_description= description,
      classifiers=[], 
      keywords= [],
      author= 'Daniela Sanchez, Aminael Sanchez',
      author_email= 'asanchez2@utpl.edu.ec',
      url= 'https://github.com/dani2s/ShadowCaster',
      license= 'GNU General Public License Version 3',
      packages= "shadowcaster",
      scripts=["bin/shadowcaster"],
      include_package_data=True,
      zip_safe=False,
      install_requires=['numpy',
                        'biopython',
                        'matplotlib',
                        'pandas',
                        'scipy',
                        'scikit-learn',
                        'ete3'], 
                        
      )

