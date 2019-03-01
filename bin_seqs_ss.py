#!/usr/bin/python3
from Bio import SeqIO
import argparse
import sys
import os
parser = argparse.ArgumentParser(description='This script bins sequences in a fasta file (slices the fasta file).  please specify all optional parameters. ')
parser.add_argument('--file', help='fasta file input')
parser.add_argument('--numOfFile', help='numOfFile')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
input_file=args_dict['file'] #extract the file option value from dict
numFile=args_dict['numOfFile']
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
#print (input_file)
def main():
	try:
		with open(input_file, "rU") as handle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')		
			Totalseqnum = len(SeqIO.to_dict(fasta_sequences))
			numEachFile = int(Totalseqnum) / int(numFile) # get the num for each file
			print("Total seq number: {}".format(Totalseqnum))
			handle.close()
	except:
		print("Unexpected error:", str(sys.exc_info()))

	try:
		if not os.path.exists("{}_parts".format(input_file)):
			os.makedirs("{}_parts".format(input_file))
		with open(input_file, "rU") as handle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')
			filecounter=1
			seqcounter=1
			filehandle=open("{}_parts/{}.part{}.fa".format(input_file,input_file,filecounter),'w')
			for entry in fasta_sequences:
				name,seq = entry.id, str(entry.seq)
				#print (name,seq)
				filehandle.write(">")
				filehandle.write(name)
				filehandle.write("\n")
				filehandle.write(seq)
				filehandle.write("\n")
				seqcounter+=1
				if seqcounter>=numEachFile and filecounter<int(numFile):  
					seqcounter=1#reset seq counter
					filecounter+=1 # increase file counter
					filehandle.close
					filehandle=open("{}_parts/{}.part{}.fa".format(input_file,input_file,filecounter),'w')#open new file to write
		
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)

		
if __name__ == "__main__": main()
