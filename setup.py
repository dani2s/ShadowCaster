#!/usr/bin/env python
from setuptools import setup

version = '0.9-beta'

description = """ShadoCaster is a software that combines two types of 
      information - sequence composition, phylogenetic likelihoods. """

setup(name='ShadowCaster',
      version= version,
      description= "Detects HGT events using a phylogenetic shadow",
      long_description= description,
      classifiers=[], 
      keywords= [],
      author= 'Daniela Sanchez, Aminael Sanchez',
      author_email= 'dasanchezsoto@gmail.com',
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
                        'xlsxwriter',
                        'scikit-learn',
                        'scipy'], 
                        #'sphinx-rtd-theme>=0.1.6',
                        #'Sphinx>=1.2.2'],
      )

