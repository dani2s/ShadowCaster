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
import shutil
import pandas as pd
from Bio import SeqIO



def phyloShadowing(orthoProb, alienOutput, queryFasta):
    imageOut = "histogram_alienLikelihoods.png"
    outputR = "alien_likelihoods.csv"
    
    os.system("phylo_shadow_model.R %s %s %s %s " %(orthoProb, alienOutput, imageOut, outputR))
    
    df1 = pd.read_csv('hgt_candidates.csv', index_col = 0)
    
    record_dict = SeqIO.to_dict(SeqIO.parse(queryFasta, "fasta"))
    with open('shadowcaster_predictions.fasta', 'w') as outWrite:
        for i in df1.index.tolist():
            outWrite.write('>%s\n%s\n' %(i,str(record_dict[i].seq)))
            
    shutil.rmtree('copy_proteomes')
    os.rmdir('sequences')
    os.remove('formatdb.log')
    os.remove('alien_svm.csv')
    

    
    
    
    
    