'''
Created on 4 ago. 2017

@author: daniela
'''
import os
import datetime
import sys
from distutils.dir_util import copy_tree

class BlastAtyipical():
    
    def __init__(self, origPath, listProteomes, alienCds, protFolder, blastp26= 'blastp', blast26db = 'formatdb'):
        self.listFiles = listProteomes
        self.alienFaa = alienCds
        self.database = blast26db
        self.blastProtein = blastp26
        self.path = origPath
        
        os.mkdir('copy_proteomes/')
        copy_tree(protFolder, 'copy_proteomes/')
        self.proteomes = 'copy_proteomes/'

        
    def createDatabase(self):
        for i in self.listFiles:
            q = os.path.join(self.proteomes, str(i))
            #print('%s -i %s.fasta -p T' %(self.database,q))
            os.system('%s -i %s -p T' %(self.database,q))
            
            
    def performBlastp(self):
        self.createDatabase()
        blast = 'blastout/'
        os.mkdir(blast)
        
        path2perl = os.path.join(self.path, 'generate_R_input.pl')
        os.system('perl %s %s %s %s %s' %(path2perl, self.alienFaa, blast, self.blastProtein, self.proteomes)) 
        #tuple, alien, blastFolder, blastp26, protFolder
 