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
parser.add_argument('-t', help='target gene ID file, leave blank for extracting all genes from GFF3 file')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
gff=args_dict['g'] #extract the file option value from dict
target=args_dict['t'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
gene_dict=defaultdict(dict) ## initialized a multi-dimension dict
target_gene_dict=defaultdict(dict)
## defined pattern that needs to be searched in the main program
GOterm_pattern = re.compile('Ontology_term=GO')
ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
GO_pattern = re.compile('Ontology_term=(.+?);')

skip_target_file_flag=0
#####################
##      main       ##
#####################
def main():
	global skip_target_file_flag
	global target
	try: 
		## go through target file
		if target is None:
			skip_target_file_flag=1
			print("target gene ID list file is not supplied, now extracting all genes from GFF file\n")
			target=gff #for filename creation
		else:
			with open(target, "rU") as handle:
				for line_raw in handle:
					#print(line_raw.strip())
					target_gene_dict[line_raw.strip()]=1
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()

	try: 
		## go through gff file
		with open(gff, "rU") as handle: 
			for line_raw in handle:
				if ( not (re.search(GOterm_pattern,line_raw) is None)):
					#print (line_raw)
					m=re.search(ID_pattern,line_raw)
					m2=re.search(GO_pattern, line_raw)					
					gene_ID='N/A'
					GOterms = 'N/A'
					if m and m2: 
						gene_ID=m.group(1)								
						GOterms=m2.group(1)
					else:
						print("Either gene ID or GO terms are not found in a line contains the text 'Ontology_term'!")
					#print(gene_ID)
					#print(GOterms)
					gene_ID=re.sub("-..$",'',gene_ID)
					GOtermsArray=re.split(",",GOterms)
					for GO in GOtermsArray:
						if gene_ID in gene_dict.keys():
							#print(GO)
							if GO in gene_dict[gene_ID].keys():
								gene_dict[gene_ID][GO]+=1
							else:
								gene_dict[gene_ID][GO]=1
						else:
							gene_dict[gene_ID][GO]=1
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()
	
	try:
		with open ("{}.GO".format(target),"w") as whandle:  # write GO terms to file
			#print(skip_target_file_flag)
			if skip_target_file_flag==0:
					for gene_ID in gene_dict:
						#print (gene_ID)
						if gene_ID in target_gene_dict:
							whandle.write(gene_ID)
							#whandle.write("\t")
							for GO in gene_dict[gene_ID]:
								whandle.write("\t")
								whandle.write(GO)
							whandle.write("\n")
			else:
				for gene_ID in gene_dict:
					whandle.write(gene_ID)
					#whandle.write("\t")
					for GO in gene_dict[gene_ID]:
						whandle.write("\t")
						whandle.write(GO)
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

