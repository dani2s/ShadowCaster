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

from Bio.Alphabet import  _verify_alphabet, IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
import os

class ValidateFiles():
    def __init__(self, qProt, qGenome, proteomeFolder):
        self.p = qProt
        self.g = qGenome
        self.pFolder = proteomeFolder
        self.proteomeList = None    
        
        
    def compareFasta(self):
        proteomeDict = SeqIO.to_dict(SeqIO.parse(self.p, "fasta"))
        genomeDict = SeqIO.to_dict(SeqIO.parse(self.g, "fasta"))
        lGeno = len(genomeDict)
        lProt = len(proteomeDict)
        
        if lGeno != lProt:
            raise ValueError("illegal query sequences, genome: %s CDS of proteome: %s" % (lGeno, lProt))
        
        
            
    def checkProtein(self):
        self.checkExtension()
        for i in self.proteomeList:
            for seq in SeqIO.parse(i, "fasta"):
                protSeq = str(seq.seq).translate(None, '*')                   
                seqObj = Seq(protSeq, IUPAC.extended_protein)
                my_seq = _verify_alphabet(seqObj)
                if my_seq is False:
                    raise TypeError("In file %s, sequence %s is not a protein" %(i, seq.id))
                
    def checkExtension(self):
        self.compareFasta()
        allFiles = os.listdir(self.pFolder) 
        self.proteomeList = []
        for a in allFiles:
            if a.endswith('.faa')  or a.endswith('.fasta'):
                wFile = os.path.join(self.pFolder, a)
                self.proteomeList.append(wFile)
            else:
                raise TypeError("File %s must be end in .fasta or .faa" %(a))
                
        extFile = os.path.splitext(self.p)[1] 
        if extFile.endswith('.faa')  or extFile.endswith('.fasta'):
            pass
        else:
            raise TypeError("File %s must be end in .fasta or .faa" %(self.p))
        self.proteomeList.append(self.p)
        

