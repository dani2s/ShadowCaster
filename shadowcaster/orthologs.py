'''
Created on 3 ago. 2017

@author: daniela
'''
from shutil import copy, rmtree, copyfile
import os
import re
import datetime

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet



class Orthologs():
    
    def __init__(self, orthoConfig, queryProt, protFolder):
        self.config = orthoConfig
        self.qProteome = queryProt
        self.proteomes = protFolder
        

    def performOrthoMcl(self):
        tempPath = os.path.join(os.getcwd(), "temporal/")
        seqPath = os.path.join(os.getcwd(), "sequences/") 
        workPath = "orthomcl/" 
        os.mkdir(seqPath)
        os.mkdir(workPath)
        print seqPath
        
        for i in os.listdir(self.proteomes):
            copy(self.qProteome, seqPath)
            
            pathProt = os.path.join(self.proteomes, i)
            copy(pathProt, seqPath)
            os.system("orthomcl-pipeline -i %s -o %s --nocompliant --yes -m %s " %(seqPath, tempPath, self.config))
            
            orthoResult = os.path.join(tempPath, 'pairs/pairs/orthologs.txt') 
            copy(orthoResult, workPath)
            #rename orthoMcl file
            pre, ext = os.path.splitext(i)
            os.rename('orthomcl/orthologs.txt', 'orthomcl/%s.txt' %pre)
            
            rmtree(tempPath)
            
            fileList = os.listdir(seqPath)
            for fileName in fileList:
                os.remove(seqPath + fileName)
        
            
"""
a = datetime.datetime.now()
os.chdir('pruebas')
confOrtho = '/home/daniela/orthomcl-pipeline/orthomcl.conf'
qProt = '/home/daniela/Proyectos/Escherichia_coli.faa'
protsFolder = '/home/daniela/Proyectos/proteomas'

test = Orthologs(confOrtho,qProt,protsFolder)
test.performOrthoMcl()


b = datetime.datetime.now()
c = b - a
print c.seconds
print c
"""