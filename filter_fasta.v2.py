#!/usr/bin/python3
from Bio import SeqIO
import argparse
import sys
from collections import defaultdict
import re
import linecache
parser = argparse.ArgumentParser(description='This script filters fasta sequences based on existence of a keyword.  please specify all optional parameters. ')
parser.add_argument('--fastafile', help='fasta file input')
parser.add_argument('--keyword', help='keyword')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
input_fastafile=args_dict['fastafile'] #extract the file option value from dict
keyword=args_dict['keyword']


if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
#print (input_file)

def main():

	try:
	#start filtering fasta file
		with open(input_fastafile, "rU") as handle, open("{}.filtered.fasta".format(input_fastafile),'w') as writehandle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name,seq, desc = entry.id, str(entry.seq), entry.description
				#print(name)
				if re.search(keyword,desc.strip()):
					writehandle.write(">{}\t{}\n{}\n".format(name,desc,seq))
		
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