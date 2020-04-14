#!/usr/bin/python3

from collections import defaultdict
import argparse
import sys
import linecache
import re
import os
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
parser.add_argument('--indir', help='input dir')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
indir=args_dict['indir'] #extract the file option value from dict

if len(sys.argv)<=1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
jobs_dict=defaultdict(dict) ## initialized a multi-dimension dict

## defined pattern that needs to be searched in the main program
gene_pattern = re.compile('sample (.+) done', flags=re.IGNORECASE)


#####################
##      main       ##
#####################	
def main():
	try: 
		###go through a dir
		with open(f'{indir}_parsed.tab', "w") as writehandle:
			writehandle.write(f"JobName\tJobID\tJobArrayNum\tSampleNum\tStatus\n")
			for filename in os.listdir(indir):
				print(f"{filename}")
				filenamefields = filename.split("_")
				with open(f"{indir}/{filename}", "r", encoding="utf-8") as handle:
					for line_raw in handle:
						m=re.search(gene_pattern,line_raw)
						if m:
							sampleNum=m.group(1)
							writehandle.write(f"{filenamefields[0]}\t{filenamefields[1]}\t{filenamefields[2].split('.')[0]}\t{sampleNum}\tDone\n")

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

