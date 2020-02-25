#!/usr/bin/python3
#this script extract GO terms associated with each AGAP identifier
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
parser= MyParser(description='this script extract GO terms associated with each AGAP identifier')
#parser = argparse.ArgumentParser()
parser.add_argument('-g', help='GFF3 file')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
gff=args_dict['g'] #extract the file option value from dict

if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
gene_dict=defaultdict(dict) ## initialized a multi-dimension dict

## defined pattern that needs to be searched in the main program
ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
desc_pattern = re.compile('description=(.+)')
gene_pattern = re.compile('\tgene\t')

#####################
##      main       ##
#####################
def main():

	try: 
		## go through gff file
		with open(gff, "rU") as handle: 
			for line_raw in handle:
				if ( not (re.search(gene_pattern,line_raw) is None)):
					#print (line_raw)
					m=re.search(ID_pattern,line_raw)
					m2=re.search(desc_pattern, line_raw)
					gene_ID='N/A'
					desc = 'N/A'
					if m: 
						gene_ID=m.group(1)								
					if m2: 
						desc=m2.group(1)

					#print(gene_ID)
					#print(desc)
					gene_ID=re.sub("-..$",'',gene_ID)
					gene_dict[gene_ID]=desc
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()
	
	try:
		with open ("{}.desc".format(gff),"w") as whandle:  # write GO terms to file

			for gene_ID in gene_dict:
				whandle.write(gene_ID)
				whandle.write("\t")
				whandle.write(gene_dict[gene_ID])
				whandle.write("\n")			
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

