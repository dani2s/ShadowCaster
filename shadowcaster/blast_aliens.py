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
 