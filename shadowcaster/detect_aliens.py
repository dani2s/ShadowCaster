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
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
    
import pandas as pd
from sklearn.svm import OneClassSVM
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import translate


def alienSVM(kmerFile, fastaCds, nuVal):
        
    #merge kmer and codon usage dataframes
    df1 = pd.read_csv('kl_codonUsage.csv', index_col = 0)
    df2 = pd.read_csv(kmerFile, index_col = 0)
    result = pd.concat([df1, df2], axis=1)
    result.columns = ['codon', 'kmer']
    
    clf = OneClassSVM(nu=float(nuVal))
    clf.fit(result)
    y_pred = clf.predict(result)
    result['pred'] = y_pred
    
    dfOutliers = result[result['pred'] == -1]
    hostDf = result[result['pred'] == 1]
    
    #Plot outliers obtained with One Class SVM
    plotSVM(hostDf, dfOutliers)
    
    #save alien sequences (protein alphabet)
    record_dict = SeqIO.to_dict(SeqIO.parse(fastaCds, "fasta"))
    with open('alien_genes.fasta', 'w') as outWrite:
        for i in dfOutliers.index.tolist():
            line = translate(str(record_dict[i].seq), to_stop=True) 
            outWrite.write('>%s\n%s\n' %(i,line))


def plotSVM(hostDf, dfOutliers):
    fig = plt.figure()
    fig.suptitle("One Class Support Vector Machine", fontsize=14)      
    ax1 = fig.add_subplot(111) #411 means 4 rows and 1 column of plot windows; the third digit specifies the order
    ax1.scatter(list(hostDf.kmer), list(hostDf.codon), label='Native genes', color = 'lavender')
    ax1.scatter(list(dfOutliers.kmer), list(dfOutliers.codon), label='Alien genes', color = 'green')
    ax1.set_ylabel('4-mers Chi-square')
    ax1.set_xlabel('Codon usage Kullback-Leibler')
    ax1.legend(loc='best', fontsize='small')
    plt.savefig('plot_aliens.png')
    #plt.show()