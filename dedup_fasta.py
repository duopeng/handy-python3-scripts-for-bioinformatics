#!/usr/bin/python3
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio import SeqIO

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
parser.add_argument('--fasta', help='fasta to dedup')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
fastafile=args_dict['fasta'] #extract the file option value from dict

if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
fasta_seqs=defaultdict(dict) ## initialized a multi-dimension dict

#####################
##      main       ##
#####################	
def main():
	try: 
		## go through fasta file
		with open(fastafile, "rU") as handle:
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name, desc,seq = entry.id, entry.description, str(entry.seq)
				fasta_seqs[seq]=name

		with  open("{}.dedup.fasta".format(fastafile.replace("\.fa","").replace("\.fasta","")),'w') as writehandle:
			for seqs in fasta_seqs.keys():
				name=fasta_seqs[seqs]
				writehandle.write(">{}\n{}\n".format(name,seqs))

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

