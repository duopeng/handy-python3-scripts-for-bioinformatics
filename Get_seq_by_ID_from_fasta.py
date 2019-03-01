#!/usr/bin/python3
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio import SeqIO

##############
## arguments##
##############
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='This script fetches seqs by ID from a fasta file.  please specify all optional parameters. ')
#parser = argparse.ArgumentParser()
parser.add_argument('--fafile', help='fasta file')
parser.add_argument('--IDfile', help='ID file,')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
fafile=args_dict['fafile'] #extract the file option value from dict
IDfile=args_dict['IDfile'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
ID_dict=defaultdict(dict) ## initialized a multi-dimension dict

#####################
##      main       ##
#####################	
def main():
	try: 
		## go through a text file
		with open(IDfile, "rU") as handle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			for line_raw in handle:
				ID_dict[line_raw.strip()]=1
		## go through fasta file
		with open(fafile, "rU") as handle, open("{}.{}.fa".format(fafile,IDfile),"w") as wh:
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name, desc,seq = entry.id, entry.description, str(entry.seq)
				flag=0
				for ID in ID_dict:
					regex_pattern = r"\b(?=\w)" + re.escape(ID) + r"\b(?!\w)" #build regex pattern
					if not (re.search(regex_pattern,name) is None):
						flag=1
				if flag==1:
					wh.write(">{}\t{}\n{}\n".format(name,desc,seq))
	
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







