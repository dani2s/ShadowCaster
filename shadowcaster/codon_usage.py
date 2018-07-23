'''
Created on 28 jun. 2017

@author: daniela
'''
from Bio import SeqIO
import pandas as pd
from collections import OrderedDict
from scipy.stats import entropy
import os
 
class CodonUsage():
    
    """Count codons and use KL as metric to find differences between all sequences and each CDS"""
    
    CodonsDict = {'TTT': 1, 'TTC': 1, 'TTA': 1, 'TTG': 1, 'CTT': 1,
                  'CTC': 1, 'CTA': 1, 'CTG': 1, 'ATT': 1, 'ATC': 1, 
                  'ATA': 1, 'ATG': 1, 'GTT': 1, 'GTC': 1, 'GTA': 1, 
                  'GTG': 1, 'TAT': 1, 'TAC': 1, 'TAA': 1, 'TAG': 1, 
                  'CAT': 1, 'CAC': 1, 'CAA': 1, 'CAG': 1, 'AAT': 1, 
                  'AAC': 1, 'AAA': 1, 'AAG': 1, 'GAT': 1, 'GAC': 1, 
                  'GAA': 1, 'GAG': 1, 'TCT': 1, 'TCC': 1, 'TCA': 1, 
                  'TCG': 1, 'CCT': 1, 'CCC': 1, 'CCA': 1, 'CCG': 1, 
                  'ACT': 1, 'ACC': 1, 'ACA': 1, 'ACG': 1, 'GCT': 1, 
                  'GCC': 1, 'GCA': 1, 'GCG': 1, 'TGT': 1, 'TGC': 1, 
                  'TGA': 1, 'TGG': 1, 'CGT': 1, 'CGC': 1, 'CGA': 1, 
                  'CGG': 1, 'AGT': 1, 'AGC': 1, 'AGA': 1, 'AGG': 1, 
                  'GGT': 1, 'GGC': 1, 'GGA': 1, 'GGG': 1}  

    SynonymousCodons = OrderedDict()
    SynonymousCodons['ALA'] = ['GCA', 'GCC', 'GCG', 'GCT']
    SynonymousCodons['CYS'] = ['TGT', 'TGC'] 
    SynonymousCodons['ASP'] = ['GAT', 'GAC']      
    SynonymousCodons['GLU'] = ['GAG', 'GAA']
    SynonymousCodons['PHE'] = ['TTT', 'TTC']
    SynonymousCodons['GLY'] = ['GGT', 'GGG', 'GGA', 'GGC']
    SynonymousCodons['HIS'] = ['CAT', 'CAC']
    SynonymousCodons['ILE'] = ['ATC', 'ATA', 'ATT']
    SynonymousCodons['LYS'] = ['AAG', 'AAA']
    SynonymousCodons['LEU'] = ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA']
    SynonymousCodons['MET'] = ['ATG'] 
    SynonymousCodons['ASN'] = ['AAC', 'AAT']
    SynonymousCodons['PRO'] = ['CCT', 'CCG', 'CCA', 'CCC']
    SynonymousCodons['GLN'] = ['CAA', 'CAG']
    SynonymousCodons['ARG'] = ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA']
    SynonymousCodons['SER'] = ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT']
    SynonymousCodons['THR'] = ['ACC', 'ACA', 'ACG', 'ACT']
    SynonymousCodons['VAL'] = ['GTA', 'GTC', 'GTG', 'GTT'] 
    SynonymousCodons['TRP'] = ['TGG'] 
    SynonymousCodons['TYR'] = ['TAT', 'TAC']
    SynonymousCodons['STOP'] = ['TAG', 'TGA', 'TAA']

    
    def __init__(self, fastaCDS):
        self.cds = fastaCDS
        self.outCds = 'codon_cds.csv'
        self.outSequence = 'codon_mergeSeq.csv'
        self.resultFile = 'kl_codonUsage.csv'
    
    def codonCount(self, seq, name):
        mergeGenome = ''
        codon_count = self.CodonsDict.copy()
        freq = OrderedDict()
      
        if str(seq).islower():
            dna_sequence = str(seq).upper()
        else:
            dna_sequence = str(seq)
           
        for i in range(0, len(dna_sequence), 3): 
            codon = dna_sequence[i:i + 3] 
            if codon in codon_count:
                mergeGenome = mergeGenome + codon 
                codon_count[codon] += 1 
            #else: 
            #    raise TypeError("illegal codon %s in gene: %s" % (codon, name)) 
        totalCodons = sum(codon_count.values())
        #print totalCodons
       
        for aa in self.SynonymousCodons: 
            total = 0.0   
            codons = self.SynonymousCodons[aa] 
            nSyn = len(codons)
         
            for codon in codons:
                total += codon_count[codon]
            freqAA = float(total)/ totalCodons 
         
            # calculate the frequencies of each codon
            for codon in codons:
                codonCount = codon_count[codon]+1
                freq[codon] = (float(codonCount)/totalCodons)*(float(nSyn)/freqAA)
      
        return freq, mergeGenome  
          


    def call_codonCount(self):
        codonsFile = OrderedDict()
        allValidecodons = ''

        #Count codons of each CDS
        for seq in SeqIO.parse(self.cds, "fasta"):
            name, seq = seq.id, str(seq.seq)
            #print name
            result, mergeLine = self.codonCount(seq, name)
            codonsFile[name] = result.values()
            allValidecodons = allValidecodons + mergeLine
      
        codonsData = pd.DataFrame.from_dict(codonsFile, orient='index', dtype=float)
        codonsData.columns = result.keys()   
        codonsData.to_csv(self.outCds) 

        #Count all codons of the complete merge sequence
        codonSeq, x = self.codonCount(allValidecodons, 'sequence')
        
        #Save merge CDS
        with open('allValidecodons.fasta', 'w') as mFile:
            mFile.write('>allValidecodons\n%s\n' %allValidecodons) 
            
        dfSeq = pd.DataFrame.from_dict(codonSeq, orient='index', dtype=float)
        newDf = dfSeq.T
        newDf.to_csv(self.outSequence)


    def metrics(self):
        self.call_codonCount()
        resultDict = OrderedDict()
        codCds = pd.read_csv(self.outCds, index_col=0)
        codSeq = pd.read_csv(self.outSequence, index_col=0)
        
        for i in codCds.index.tolist():
            
            #var1 corresponds to frequencies from genes and var2 from the complete sequence 
            var1 = codCds.loc[i,:].values
            var2 = codSeq.iloc[0,:].values
            kl = entropy(var1, var2)
            resultDict[i] = kl
          
        scoresDf = pd.DataFrame.from_dict(resultDict, orient='index', dtype=float)
        scoresDf.to_csv(self.resultFile)
        os.remove('allValidecodons.fasta')
        
