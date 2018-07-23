'''
Created on 10 oct. 2017

@author: daniela
'''
import os
import datetime

def phyloShadowing(orgPath, orthoProb, alienOutput):
    imageOut = "histogram_likelihood.png"
    outputR = "likelihoods.csv"
    
    path2Rscript = os.path.join(orgPath, 'R/phylo_shadow_model.R')
    os.system("Rscript %s %s %s %s %s " %(path2Rscript, orthoProb, alienOutput, imageOut, outputR))
    
"""
a = datetime.datetime.now()
orgPAth = '/home/daniela/workspace/reports/'
aFil = "/home/daniela/Proyectos/Hemme_2016/orthologs_probabilities15.csv"
bFile = "/home/daniela/workspace/reports/2018-01-15_16-41-57/phylogenetic/R_alien.csv"
phyloShadowing(orgPAth,aFil, bFile)

b = datetime.datetime.now()
c = b - a
print c.seconds
print c
"""