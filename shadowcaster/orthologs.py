'''
    This file is part of ShadowCaster.
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
        
            