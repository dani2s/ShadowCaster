'''
Created on 10 oct. 2017

@author: daniela
'''
from Bio import SeqIO
from argparse import ArgumentParser
import ConfigParser


def arguments():
    parser = ArgumentParser()
    config = ConfigParser.ConfigParser()
    
    parser.add_argument('--config_file') 
    args = parser.parse_args()
    
    #Read config_file
    config.read(args.config_file)
    
    #Create a dictionary to store values of the config_file
    config_d = {}
    #COMPLETE PATH

    config_d['qGenome'] = str(config.get('Files', 'query_genome'))
    config_d['qProteome'] = str(config.get('Files', 'query_proteome'))
    config_d['ortho'] = str(config.get('Files', 'orthomcl_config'))
    
    config_d['prots'] = str(config.get('Path', 'proteomes_folder'))
    config_d['blastp'] = str(config.get('Path', 'blastp26'))
    config_d['db26'] = str(config.get('Path', 'formatdb26'))

    #config_d['kmer_len'] = str(config.get('Parametric', 'kmer_length'))
    #config_d['metric'] = str(config.get('Parametric', 'metrics'))

    config_d['nuSVM'] = str(config.get('Phylogenetic', 'quantile'))
    
    return config_d



