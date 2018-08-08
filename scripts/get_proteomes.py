#!/usr/bin/env python
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
from ete3 import NCBITaxa
from argparse import ArgumentParser
import pandas as pd

import os
import sys
import datetime
import subprocess
import glob

def arguments():
    parser = ArgumentParser()

    #name of the query specie
    parser.add_argument('-n', '--sp_name',
        help=("Specify name of the query sp, just two words.\
        Between the two words there MUST be an underscore(Escherichia_coli)"))
    
    #Handle number of species for the shadow
    parser.add_argument('-sp', '--num_sp', default=25, type=int,
      help='Specify number of species to build the phylogenetic shadow, default 25.')
    
    #update taxonomy database
    parser.add_argument('-db', '--update_db', default=False, action='store_true',
      help=(" First time using? ETE package will download the latest NCBI taxonomy database.\
      If True, ETE will download the latest db from NCBI and overwrite the current local database, default False"))      
    #Output path
    parser.add_argument('-o', '--out_dir', 
      help=("Specify an existing directory where output will be placed."
            "If not provided, a directory with the sp name provided will be created."))
    
    args  = parser.parse_args()
        
    return args


def check():
    rc = subprocess.call(['which', 'esearch'])
    if rc == 0:
        print 'esearch installed!\n'
    else:
        print 'esearch missing in path! Please check the path in .profile or .bash_profile'
        sys.exit(0)


def get_taxonomy(updateBool, spName):
    ncbi = NCBITaxa()
    
    #add update condition
    if updateBool is True:
        ncbi.update_taxonomy_database()
    
    #get only genus name
    genus = spName.partition('_')[0]
    
    name2taxid = ncbi.get_name_translator([genus])
    
    lineage = ncbi.get_lineage(name2taxid[genus][0])
    
    return lineage[2:]

def get_ids(url, exclude_name, maxNum, ids4download):
    
    df1 = pd.read_table(url, sep='\t', header=None)
    df1.columns = ['name', 'url']
    #print list(df1.name)
    for i in df1.index.tolist():
        lenDict = len(ids4download)
        valName = df1.loc[i,'name']
        valUrl = df1.loc[i,'url']
        if lenDict < maxNum:
            nameComplete = valName.split()
            binomialName = nameComplete[0] + '_' + nameComplete[1]
            #print binomialName
            if binomialName != exclude_name and binomialName not in ids4download.keys():
                ids4download[binomialName] = valUrl
        else:
            break
         
    return ids4download

def get_ftpPath(dictSp):
    ncbiComnd= '''esearch -db assembly -query '%s' | 
    efetch -format docsum | xtract -pattern DocumentSummary -element Organism -block FtpPath -match "@type:RefSeq" -element FtpPath | 
    sed 's/$/\//' > %s'''
    
    #Get ftp paths for family and order of query sp(include reference and representative genomes)
    for i in ['family', 'order', 'otherkingdom']:
        txid = dictSp[i][2]
        searchArg1 = '''txid%s[organism] AND "complete genome"[filter] AND "latest refseq"[Filter] AND ("reference genome"[Filter] OR "representative genome"[Filter])''' %txid
        os.system(ncbiComnd %(searchArg1, dictSp[i][1]))
    
    #Get ftp paths for phylum and kingdom of query sp (only reference genomes)  
    for j in ['phylum', 'kingdom']:
        txid1 = dictSp[j][2]
        searchArg2 = '''txid%s[organism] AND "complete genome"[filter] AND "latest refseq"[Filter] AND "reference genome"[Filter]''' %txid1
        os.system(ncbiComnd %(searchArg2, dictSp[j][1]))

        
def get_ftpPath4archaea(dictSp):
    #Only for archae, different search term for taxonomy ranks (include reference and representative genomes)
    
    ncbiComnd= '''esearch -db assembly -query '%s' | 
    efetch -format docsum | xtract -pattern DocumentSummary -element Organism -block FtpPath -match "@type:RefSeq" -element FtpPath | 
    sed 's/$/\//' > %s'''
    
    #Get ftp paths 
    for a in ['family', 'order', 'phylum', 'kingdom']:
        txid = dictSp[a][2]
        searchArg1 = '''txid%s[organism] AND "complete genome"[filter] AND "latest refseq"[Filter] AND ("reference genome"[Filter] OR "representative genome"[Filter])''' %txid
        os.system(ncbiComnd %(searchArg1, dictSp[a][1]))
    
    #Get ftp paths for other superkingdom of query sp    
    for b in ['otherkingdom']:
        txid1 = dictSp[b][2]
        searchArg2 = '''txid%s[organism] AND "complete genome"[filter] AND "latest refseq"[Filter] AND "reference genome"[Filter]''' %txid1
        os.system(ncbiComnd %(searchArg2, dictSp[b][1]))
    
    
def main(args):
    if args.out_dir == None: 
        os.mkdir('%s/' %args.sp_name)
        os.chdir('%s/' %args.sp_name)
    else:
        os.mkdir('%s' %args.out_dir)
        os.chdir('%s' %args.out_dir)
        
          
    lineage = get_taxonomy(args.update_db, args.sp_name)
    #print lineage
    
    numIds = int(args.num_sp * 0.20)
    otherKngdm = (4*numIds) + int(args.num_sp * 0.08)
        
    
    if lineage[0] == 2:
        txOtherKingdm = 2157
    else:
        txOtherKingdm = 2
    
    mainDict = {'family': [2*numIds, 'ftpdirpaths_family.txt', lineage[4]],
                'order': [3*numIds, 'ftpdirpaths_order.txt', lineage[3]],
                'phylum': [4*numIds, 'ftpdirpaths_phylum.txt', lineage[1]],
                'kingdom': [int(args.num_sp), 'ftpdirpaths_reino.txt', lineage[0]],
                'otherkingdom': [otherKngdm, 'ftpdirpaths_otroreino.txt', txOtherKingdm]}
    
    orden = ('family', 'order', 'phylum', 'otherkingdom', 'kingdom')
    
    #Get ftp paths of each rank to download
    if mainDict['kingdom'][2] == 2:
        get_ftpPath(mainDict)
    else:
        get_ftpPath4archaea(mainDict)
    
    #Filter organisms according to the weight of the shadow
    idsFinal = {}
    
    for i in orden:
        idsFinal = get_ids(mainDict[i][1], args.sp_name, mainDict[i][0], idsFinal)
    
    #Write log fie
    with open('log.txt', 'w') as outFile:
        outFile.write('Downloaded species can be found in these ftp urls\n')
        for it in idsFinal.keys():
            outFile.write('%s\t%s\n' %(it,idsFinal[it]))
    
    #remove ftp address searched
    for akeys in mainDict.keys():
        os.remove(mainDict[akeys][1])
        
    os.mkdir('proteomes/')
    print "Downloading from FTP... \nMay take several minutes \n"
    for p in idsFinal.keys():
        url = idsFinal[p]
        #print p, url
        os.system('wget -q -c %s*_protein.faa.gz' %url)
        
        fileFaa = glob.glob('*_protein.faa.gz')
        if os.path.isfile(fileFaa[0]):
            os.system('gzip -d *_protein.faa.gz')
            os.system('mv *_protein.faa proteomes/%s.fasta' %p)
        else:
            print 'Cannot connect with NCBI ftp server\n'
            print 'The log file contains the name of the species and its ftp address\n'
            sys.exit(0)
    
    print 'The log file contains the name of the downloaded species and its ftp address\n'
    print 'The proteomes folder contains fasta files to run ShadowCaster, please copy the full path in the configuration file\n'
            
        
if __name__ == '__main__':
    args = arguments()
    check()
    a = datetime.datetime.now()
    main(args)
    b = datetime.datetime.now()
    c = b - a
    print 'Took %s to complete' %c   
