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

from argparse import ArgumentParser
import ConfigParser

def get_version():
  from shadowcaster import __version__
  return '%(prog)s {version}'.format(version=__version__)

def arguments():
    parser = ArgumentParser()
    config = ConfigParser.ConfigParser()
    
    parser.add_argument('--config_file', 
        help=("Specify the configuration file, containing all the paths \
            needed to run ShadowCaster. A template of these file(args.ini) can be find \
            in the bin folder."))

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

    config_d['nuSVM'] = float(config.get('Phylogenetic', 'quantile'))

    parser.add_argument('-v','--version', action='version',
      version=get_version())
    
    return config_d



