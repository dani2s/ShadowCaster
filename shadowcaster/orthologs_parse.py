'''
Created on 11 jul. 2017

@author: daniela
'''
import re
import os
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
import pandas as pd
import numpy as np

import datetime

class ParseOrthologsResults():
    
    def __init__(self, proteomesFolder, queryProteome):
        self.workPath = 'orthomcl/'
        self.protFolder = proteomesFolder
        self.queryProt = queryProteome
        self.listFiles = None
        
        
    def getFiles(self):
        faa_files = [f for f in os.listdir(self.protFolder) if f.endswith('.faa')]
        fasta_files = [f for f in os.listdir(self.protFolder) if f.endswith('.fasta')]
        self.listFiles = fasta_files + faa_files               
        
        
    def getIdentity(self):
        self.getFiles()
        needleOut = self.workPath + 'needle.txt'
        dfResult = 'orthologs_probabilities.csv'
                
        patt = re.compile(r'.*\|(.*)\s+.*\|(.*)\s+\d+')
        pattFile = re.compile(r'\.')
        identPatt = re.compile(r'.*Identity\:\s+(\d+)\/(\d+)')
        
        hostSeq = SeqIO.to_dict(SeqIO.parse(self.queryProt, "fasta"))
        hostCds = len(hostSeq)
        #print hostCds
        
        with open(dfResult, 'w') as handle:   
            handle.write('name\torthologs\tmean\tstd\tPr_orthologs\n')
            
            for i in range(0, len(self.listFiles)):
                
                #Define variable
                qSeq = os.path.join(self.protFolder, self.listFiles[i])
                queryGen = SeqIO.to_dict(SeqIO.parse(qSeq, "fasta"))
                name = str(pattFile.split(self.listFiles[i])[0])
                pairsFile = self.workPath + name + '.txt'
                
                queryID = []
                hostID = []
                ident = []
                
                for line in open(pairsFile, 'r'):
                    objM = patt.match(line)
                   
                    if objM.group(1)in hostSeq.keys():
                        hostProt = objM.group(1)
                        qProt = objM.group(2)
                    else:
                        hostProt = objM.group(2)
                        qProt = objM.group(1)
                        
                    needle_cline = NeedleCommandline(asequence= 'asis:%s' %str(hostSeq[hostProt].seq), 
                                                    bsequence= 'asis:%s' %str(queryGen[qProt].seq), 
                                                    gapopen=10, gapextend=0.5, outfile= needleOut)
                    needle_cline()
                        
                        
                    for liFile in open(needleOut, 'r'):
                        objIdent = identPatt.match(liFile)
                        if objIdent is not None:
                            scoreIdent = float(objIdent.group(1))/float(objIdent.group(2))
                            #print scoreIdent
                            hostID.append(hostProt)
                            queryID.append(qProt)
                            ident.append(scoreIdent)
                                
                #print len(ident), len(coliID), len(queryID)
                #Identities of each organism
                #ortList = [('queryProteome_ID',hostID),('%s_ID' %name,queryID),('Identity',ident)] 
                orthoTotal = len(ident)
                #print orthoTotal
                #orthoData = pd.DataFrame.from_items(ortList)
                #orthoData.to_csv('ortho_%s.csv' %name, index=False)
                a = np.array(ident)
                meanQuery = np.mean(a)
                stdQuery = np.std(a)
                #print hostCds
                orthoPr = float(orthoTotal)/hostCds
                handle.write('%s\t%s\t%s\t%s\t%f\n' %(name, orthoTotal, meanQuery, stdQuery, orthoPr))    
                #print '%s\t%s\t%s\t%s\t%f\n' %(name, orthoTotal, meanQuery, stdQuery, orthoPr)
        os.remove(needleOut)
        return self.listFiles

