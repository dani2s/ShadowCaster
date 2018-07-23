#!/usr/bin/env python
'''
Created on 21 may. 2018

@author: Daniela Sanchez
'''

from ete3 import NCBITaxa
from Bio import Entrez, SeqIO
from argparse import ArgumentParser

import os
import datetime
import re
import random


def arguments():
    parser = ArgumentParser()

    #Input file
    parser.add_argument('proteome_file', 
        help=("Specify the proteome file of the query specie, containing sequences in fasta format."))
    
    #name of the query specie
    parser.add_argument('sp_name',
        help=("Specify name of the query specie, just two words.\
        Between the two words there MUST be an underscore(Escherichia_coli)"))
    
    #email for NCBI
    parser.add_argument('-e','--email',
        help='Email address to fetch data through NCBI')
    
    #Handle number of species for the shadow
    parser.add_argument('-sp', '--num_sp', default=25, type=int,
      help='Specify number of species to build the phylogenetic shadow, default 25.')
    
    #update taxonomy database
    parser.add_argument('-db', '--update_db', default=False, action='store_true',
      help=(" First time using? ETE package will download the latest NCBI taxonomy database.\
      If True, ETE will download the latest db from NCBI and overwrite the current local database, default False"))      
    #Output path
    parser.add_argument('-o', '--out_dir', default=os.curdir,
      help=("Specify an existing directory where output will be placed."
            "If not provided, current directory will be used."))
    
    args  = parser.parse_args()
        
    return args
    


#check the input check check



def get_taxonomy(updateBool, spName):
    ncbi = NCBITaxa()
    
    #add update condition
    if updateBool is True:
        ncbi.update_taxonomy_database()
    
    #get only genus name
    genus = spName.partition('_')[0]
    
    name2taxid = ncbi.get_name_translator([genus])
    #print name2taxid[genus][0]
    lineage = ncbi.get_lineage(name2taxid[genus][0])
    
    return lineage[2:]

def ncbiSearch(mail,dbase, termSearch, retMax):
    Entrez.email = mail
    handle = Entrez.esearch(db=dbase, term = termSearch, retmax=retMax)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def getIds(mail,lineage, spName, numGenus, numRank, numKing, numKingOther):  
        
    name4search = spName.replace("_", " ") 
    
    #Get Genus organisms
    #numberGenus= numGenus + 1  #add 1 to exclude the same sp later
    genusSearch = ncbiSearch(mail, "nucleotide", "(txid%s[Organism] AND complete genome[title]" %lineage[-1], numGenus)
    excludeOrg = ncbiSearch(mail, "nucleotide", "%s[Organism] AND complete genome[title]" %str(name4search), numGenus)
    genusOrg = list(set(genusSearch) - set(excludeOrg))
    
        
    #Get ids for superkingdoms
    superKingdomOrg = ncbiSearch(mail, "nucleotide", "(txid%s[Organism] AND complete genome[title]" %lineage[0], numKing)
    #Get ids for archae
    if lineage[0] == 2:
        superKingdomOther = ncbiSearch(mail, "nucleotide", "(txid2157[Organism] AND complete genome[title]", numKingOther)
    #Get ids for bacteria
    if lineage[0] == 2157:
        superKingdomOther = ncbiSearch(mail, "nucleotide", "(txid2[Organism] AND complete genome[title]", numKingOther)
    
    
    supKg = superKingdomOrg + superKingdomOther
    
    
    #Get sequences for other taxonomic ranks
    rankOrgs = []
    for i in lineage[1:len(lineage)-1]:
        termSearch = "(txid%s[Organism] AND complete genome[title]" %i
        resultsNcbi = ncbiSearch(mail, "nucleotide", termSearch, numRank*2)
        #print i , len(resultsNcbi)
        for val in resultsNcbi:
            if val not in rankOrgs:
                rankOrgs.append(val)       
    
    return excludeOrg, genusOrg, rankOrgs, supKg



def ncbiFetch(mail, idSp):
    Entrez.email = mail
    #Get the organism name to save as org_name.faa
    handle = Entrez.esummary(db="nucleotide", id=str(idSp))
    record = Entrez.read(handle)
    nameRecord = str(record[0]['Title']).split()
    faaName = nameRecord[0] + '_' + nameRecord[1] + '_' + idSp
    print str(record[0]['Title']),'\tid:', idSp 
    #faaName = nameRecord.replace(" ", "_") 
    handle.close()
    
    #Fetch all the cds(aa) of the specie
    handle1 = Entrez.efetch(db="nucleotide", id =str(idSp), rettype="fasta_cds_aa", retmode="text")
        
    #Save each sequence only with protein id(OrthoMCl requirements)
    pattern = re.compile(r'.*\[protein_id=(\w+\d+.\d)].*')
    ids = []
    with open("proteomes/%s.fasta" %faaName, 'w') as fafile:
        for seq_record in SeqIO.parse(handle1, "fasta"):
            pattResult = pattern.search(seq_record.description)
            if not pattResult is None:
                protId = str(pattResult.group(1))
                #Write only unique ids
                if protId not in ids:
                    fafile.write(">%s\n%s\n" %(protId,str(seq_record.seq)))
                    ids.append(protId)
    handle1.close()
        
  
def main(args):    
    os.mkdir('proteomes/')
    lineage = get_taxonomy(args.update_db, args.sp_name)
    #print lineage
    
    #define numbers of sequences for the search
    blast = int(args.num_sp * 0.20)
    genus = int(args.num_sp * 0.20)
    rankTax = int(args.num_sp * 0.40)
    ownSuperkingdom = int(args.num_sp * 0.12)
    otherSuperkingdom = int(args.num_sp * 0.08)
    
    identOrg, genusList, rankList, kingList = getIds(args.email,lineage, args.sp_name, genus, rankTax, ownSuperkingdom, otherSuperkingdom)
    
    #blast results
    #blastList
     
    #Remove duplicate from the other lists
    list1 = genusList + kingList
    lenList1 = len(list1) + len(identOrg)
    rank = list(set(rankList) - set(list1))   
    
    #rankTax plus the missing ids of other lists
    sampleRank = int(args.num_sp) - lenList1
    rankIds = random.sample(rank, sampleRank)
    finalIds = rankIds + list1
    
    #print finalIds
    print 'Downloading cds(aa) of these organisms:'
    for i in finalIds:
        ncbiFetch(args.email, i)
    print "Name of each fasta file: organism's binomial name and id(nucleotide database)"
    
    #print "END :)"
    
if __name__ == '__main__':
    args = arguments()
    a = datetime.datetime.now()
    main(args)
    b = datetime.datetime.now()
    c = b - a
    print 'Took %s to complete' %c
