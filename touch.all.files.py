#touches all files (including subdirectories) of PATH
#usage:  python touch.all.files.py --path [PATH]
import os
import glob
import argparse
import sys
import linecache
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
parser.add_argument('--path', help='path')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
path=''
path=args_dict['path'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
	path = '/n/scratch2/dp235'
	print("no path argument supplied, using default path /n/scratch2/dp235")
#    parser.print_help()
#    sys.exit(1)

def main():
	try: 

		##list files in path

		for p, subdirs, files in os.walk(path):
			for name in files:
				file_path=os.path.join(p, name)
				print(file_path)
				subprocess.call("touch {}".format(file_path), shell=True ) # wait for command to finish
			
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
