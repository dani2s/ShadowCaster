#!/usr/bin/env python
'''
    ShadowCaster
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2018  Daniela Sanchez and Aminael Sanchez

    ShadowCaster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ShadowCaster is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Questions? Please contact asanchez2utpl.edu.ec

'''
import os
import datetime
import sys
import inspect

from shadowcaster.parser import arguments
from shadowcaster.parametric import ParametricMethods
from shadowcaster.codon_usage import CodonUsage
from shadowcaster.detect_aliens import alienSVM
from shadowcaster.orthologs import Orthologs
from shadowcaster.orthologs_parse import ParseOrthologsResults
from shadowcaster.blast_aliens import BlastAtyipical
from shadowcaster.likelihood import phyloShadowing
import shadowcaster.likelihood


def main(config):
    
    
    originalPath= os.path.dirname(inspect.getfile(shadowcaster.likelihood))


    #Make working directory
    wDir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    os.mkdir(wDir)      
    os.chdir(wDir)
        
    #Create log file
    logFile = open('log.txt', 'w')
    logFile.write('Arguments of ShadowCaster \nConfiguration file:\n')
    for key, value in config.iteritems():
        logFile.write('%s %s\n' %(key,value))
    logFile.close()
        
    
    print >> sys.stderr, "Step 1: Checking input files\n"
       

    #Create folder for parametric's results
    os.mkdir('parametric/')
    os.chdir('parametric/')
    phase1 = os.getcwd()

    print >> sys.stderr, "Step 2: Starting Composition Methods"
    
    #Instance Parametric class  (13sec)
    faaMerge = 'mergeCDS.fasta' 
    kmer = ParametricMethods(config['qGenome'], faaMerge,'chi2',4,reverse=False)
    kmerDf = kmer.metricMeasure()
    print >> sys.stderr, "=====k-mer frequecies are calculated===="

    #Instance Codon Usage class (5 sec)
    CU = CodonUsage(config['qGenome'])
    CU.metrics()
    print >> sys.stderr, "=======Codon Usage is calculated======="
    
    
    #Call R script  (6 sec)
    print >> sys.stderr, "===Detecting alien with One Class SVM==="
    phase1Result = '%s/alien_genes.fasta' %phase1
    alienSVM(kmerDf, config['qGenome'], config['nuSVM'])


    print >> sys.stderr, "==Composition results can be found in parametric folder==\n"
 
    os.chdir('..')
    
    
    os.mkdir('phylogenetic')
    os.chdir('phylogenetic')
    
    print >> sys.stderr, "Step 3: Starting Phylogenetic Methods\n"
    
    print >> sys.stderr, "====Starting OrthoMcl===="
    #Instance Orthologs class (1 hour each round 10-15 min)
    orthoMcl = Orthologs(config['ortho'], config['qProteome'], config['prots'])
    orthoMcl.performOrthoMcl() 
    print >> sys.stderr, "====OrthoMcl finished==="

    #Instance Parse orthologs results (13 min)
    orthoResults = ParseOrthologsResults(config['prots'], config['qProteome'])
    listProteomes = orthoResults.getIdentity()
    print listProteomes
    print >> sys.stderr, "====Orthologs results are parsed===\n"
    
    #Instance Blast atypical class (5 min)
    atypical = BlastAtyipical(originalPath,listProteomes, phase1Result, config['prots'], config['blastp'], config['db26'])
    atypical.performBlastp()
    print >> sys.stderr, "====Blast atypical genes done====\n\n"
    
    #Function for calculating likelihood (1min)
    #'/home/daniela/Proyectos/Hemme_2016/orthologs_30.csv'
    phyloShadowing(originalPath, 'orthologs_probabilities.csv', 'alien_svm.csv')
    
    
    print >> sys.stderr, "===Phylogenetic results done===\n\n"


if __name__ == '__main__':
    a = datetime.datetime.now()
    print "******************************"
    print "*********ShadowCaster*********"
    print '******************************\n\n'
    print 'Date:%s' %a
    config = arguments()
    results = main(config)
    b = datetime.datetime.now()
    c = b - a
    #print c.seconds
    print >> sys.stderr, "\n\nShadowCaster finished"
    print 'Took %s to complete' %c
    

