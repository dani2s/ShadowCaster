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
from collections import OrderedDict
from numpy.linalg import inv
from itertools import product, tee, izip
from Bio import SeqIO
from scipy.spatial import distance
from scipy.stats import chisquare, entropy, pearsonr
import numpy as np
import pandas as pd
import re

import datetime


class ParametricMethods():
    
    def __init__(self, windowsFile, genomeFile, metricMeasure, kmer = 4, reverse = True):
        self.windows = windowsFile
        self.genome = genomeFile
        self.k_length = kmer
        self.metric = metricMeasure
        self.revComp = reverse
        #self.folder = outPath
        self.kmerWin = None
        self.kmerGen = None

    
    def generateDictionary(self, kmer_len):
        #Generate k-mer dictionary, k is defined by the user
        baseComplement = {"A":"T","T":"A","G":"C","C":"G"}
        kmerDict = {}
        counter = 0
        for kmer in product("ATGC",repeat=kmer_len):
            if kmer not in kmerDict:
                kmerDict[kmer] = counter
                if self.revComp is True:
                    rev_compl = tuple([baseComplement[x] for x in reversed(kmer)])
                    kmerDict[rev_compl] = counter
                counter += 1
        return kmerDict, counter   
    
    def createMergeCDS(self):
        l = ''
        with open(self.genome, 'w') as mFile:
            for ss in SeqIO.parse(self.windows, 'fasta'):
                l = l + ss.seq
            mFile.write('>merge\n%s\n' %l) 
                
     
    def calculateKmers(self, kmer_hash, counter,fastaFile):
        #Calculate frequencies of each k-mer per window
        composition_d = OrderedDict()
        for seq in SeqIO.parse(fastaFile, "fasta"):
            kmers = [
                kmer_hash[kmer_tuple]
                for kmer_tuple 
                in self.window(str(seq.seq))
                if kmer_tuple in kmer_hash
                ]
            kmers.append(counter - 1)
            composition_v = np.bincount(np.array(kmers))
            composition_v[-1] -= 1
            composition_d[seq.id] = composition_v + np.ones(counter)
        composition_data = pd.DataFrame.from_dict(composition_d, orient='index', dtype=float) 
        return composition_data
                 
    def window(self, str_seq):
        #Generate substrings of k length
        els = tee(str_seq, self.k_length)
        for i,el in enumerate(els):
            for _ in xrange(i):
                next(el, None)
        return izip(*els)   
   
    def normalizeData(self, kmerData):
        #Normalize kmer frequencies, ln(p_ij) = ln[(X_ij +1) / rowSum(X_ij+1)]
        compositionDf= np.log(kmerData.divide(kmerData.sum(axis=1),axis=0))
        return compositionDf
          
    def mahalanobisDist(self,winFreq, genomeFreq):
        #Calculate the inverse of the covariance matrix
        invCov = inv(np.cov(winFreq, genomeFreq))
        print invCov[0][1] 
        #Calculate mahalanobis distance
        mahaDist = distance.mahalanobis(winFreq, genomeFreq, invCov[0][1])
        #return mahaDist

    def correlation(self,winFreq, genomeFreq):
        #Calculate Pearson correlation
        coefCor, pvalue = pearsonr(winFreq, genomeFreq)
        return coefCor

    def covariance(self,winFreq, genomeFreq):
        #Calculate covariance distance
        dimCol = len(winFreq)
        covD = np.divide(np.sum(np.dot(winFreq, genomeFreq)), dimCol)
        return covD
            
    def chiSquare(self,winFreq, genomeFreq):
        #Calculate a one-way chi square test
        chiSq , pvalue = chisquare(winFreq, genomeFreq)
        return chiSq

    def euclideanDist(self,winFreq, genomeFreq):
        #Calculate the Euclidean distance between two 1-D arrays
        euDist = distance.euclidean(winFreq, genomeFreq)
        return euDist
    
    def klDivergenge(self,winFreq, genomeFreq):
        #Calculate the Kullback-Leibler divergence 
        kl = entropy(winFreq, genomeFreq)
        return kl
    
    def delta(self,winFreq, genomeFreq):
        #Calculate the difference (delta) of a window in the complete sequence
        dimCol = len(winFreq)
        dDiff = np.divide(np.absolute(np.sum(np.subtract(winFreq, genomeFreq))), dimCol)
        return dDiff

    def generateData(self):
        #generate k-mer dataframe of the windows and genome file
        self.createMergeCDS()
        
        kmerDict, counter = self.generateDictionary(self.k_length)
        dataWin = self.calculateKmers(kmerDict, counter, self.windows)
        #For euclidean distance use raw frequencies
        if self.metric != 'euclidean':
            self.kmerWin = self.normalizeData(dataWin)
        else:
            self.kmerWin = dataWin
        self.kmerWin.to_csv('data%smerCds.csv' %(self.k_length))
        
        dataGen = self.calculateKmers(kmerDict, counter, self.genome)
        if self.metric != 'euclidean':
            self.kmerGen = self.normalizeData(dataGen)
        else:
            self.kmerGen = dataGen
        self.kmerGen.to_csv('data%smerGen.csv' %(self.k_length))
        
    
    def metricMeasure(self):
        self.generateData()
        #Measure the distance/score of window to/in complete sequence
        #pattern = re.compile(r'_')
        resultDict = OrderedDict()
        metricDict = {'mahalanobis': self.mahalanobisDist, 'correlation': self.correlation, 'covariance': self.covariance,
                      'chi2': self.chiSquare, 'euclidean': self.euclideanDist, 'Kullback-Leibler': self.klDivergenge,
                      'delta' :self.delta}
        
        for i in self.kmerWin.index.tolist():
            #index1 = str(pattern.split(i)[0])
            #var1 corresponds to frequencies from windows and var2 from the complete sequence 
            var1 = self.kmerWin.loc[i,:].values
            var2 = self.kmerGen.iloc[0,:].values
            
            resultMetric = metricDict[self.metric](var1, var2)
            resultDict[i] = resultMetric
        scoresDf = pd.DataFrame.from_dict(resultDict, orient='index', dtype=float)
        outHandle = 'data%smer_%s.csv' %(self.k_length,self.metric)
        scoresDf.to_csv(outHandle)
        return outHandle

