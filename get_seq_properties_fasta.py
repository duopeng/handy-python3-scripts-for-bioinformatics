#!/usr/bin/python3
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio.Seq import Seq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import collections
import subprocess
##############
## arguments##
##############
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='This script does xxx')
#parser = argparse.ArgumentParser()
parser.add_argument('--fasta', help='fasta file')
parser.add_argument('--len_cutoff', help='ORF length cutoff')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
fastafile=args_dict['fasta'] #extract the file option value from dict
len_cutoff=args_dict['len_cutoff']
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
list_dict=defaultdict(dict) ## initialized a multi-dimension dict
dict_genes_w_ATG
dict_genes_ORF_length


## defined pattern that needs to be searched in the main program
#ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
#parent_pattern = re.compile('Parent=(.+?)(?:\.|:|;|$)')
ATG_pattern=re.compile('^ATG')

#####################
##      main       ##
#####################	
def main():
	try: 

		## go through fasta file
		with open(fastafile, "rU") as handle, open("{}.tab".format(fastafile),'w') as writehandle: 
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name, desc,seq = entry.id, entry.description, str(entry.seq)
				#start with ATG
				if re.search(ATG_pattern,seq,flags=re.IGNORECASE):
					dict_genes_w_ATG[name]=1
				#ORF length
				coding_dna = Seq(seq, IUPAC.unambiguous_dna)
				translation = coding_dna.translate(to_stop=True)
				dict_genes_ORF_length[name]=len(translation)*3
				print(translation)
				

	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()
		
##########################
## function definitions ##
##########################
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

	
if __name__ == "__main__": main()	

