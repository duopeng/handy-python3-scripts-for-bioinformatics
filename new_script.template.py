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
parser.add_argument('--file1', help='file1')
parser.add_argument('--file2', help='file2,')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
file1=args_dict['file1'] #extract the file option value from dict
file2=args_dict['file2'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
list_dict=defaultdict(dict) ## initialized a multi-dimension dict

## defined pattern that needs to be searched in the main program
ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
parent_pattern = re.compile('Parent=(.+?)(?:\.|:|;|$)')

global myGlobalVariable

#####################
##      main       ##
#####################	
def main():
	try: 
		nonlocal myGlobalVariable
		myGlobalVariable=1 
		global myGlobalVariable
		myGlobalVariable=1
		# or else another local variable named myGlobalVariable will be created
		#To rebind variables found outside of the innermost scope, the nonlocal statement can be used; if not declared nonlocal, those variables are read-only (an attempt to write to such a variable will simply create a new local variable in the innermost scope, leaving the identically named outer variable unchanged).
		
		## go through a text file
		with open(input_file, "rU") as handle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			for line_raw in handle:
				lineSplitArray=re.split(",",line_raw.strip())
		## go through fasta file
		with open(input_fastafile, "rU") as handle, open("{}.primers.copynum.tab".format(input_fastafile),'w') as writehandle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name, desc,seq = entry.id, entry.description, str(entry.seq)
	
		##reverse complement
		from Bio.Seq import Seq
		seq=str(Seq(seq).reverse_complement())
		# GGGCCCGA
		
		##translation
		from Bio.Seq import Seq
		from Bio.Alphabet import IUPAC
		coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
		coding_dna.translate(to_stop=True)
		
		##find text occurrence
		patterntext=re.escape("text")
		patterntext=re.compile("text")
		Fprimer_occurence=len(re.findall(patterntext, string ,flags=re.IGNORECASE))
		##Or
		if (re.search("NNNNNNNNNN",seq[extract_start:extract_end+1],flags=re.IGNORECASE) is None):
		
		##substitution
		pattern.sub(replacement,string)
		
		## extract number from desc=">CurContig109_quiver_pilon_pilon_pilon_pilon_691221_694381_1 | FPrate:0.100 | OMEGA:H-1027"
		FPrate = desc.split("|")[1].lstrip().rstrip().split(":")[1]
		
		##sort dict based on value
		for match in sorted(match_dict, key=match_dict.get, reverse=True):
		
		##sort dict based on length of key
		for seq in sorted(grna_dict, key=len, reverse=True):
		
		###go through a dir
		for filename in os.listdir(input_folder):
		
		###make dir
		os.mkdir
		
		###check dir
		os.path.isdir("/home/el")
		os.path.exists("/home/el/myfile.txt")
		
		## call system command
		subprocess.call("python3 find_primer_copynum.py --fa allRepeats.fa --p {}/{} &".format(indir,filename), shell=True ) # run in background
		subprocess.call("python3 find_primer_copynum.py --fa allRepeats.fa --p {}/{}".format(indir,filename), shell=True ) # wait for command to finish		
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

