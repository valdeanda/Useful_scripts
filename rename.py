import sys
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--mapping-file", metavar="F",
      type=argparse.FileType("r"), required=True,
      help="A csv mapping file containing 2 columns - 'old-name,new-name'")
parser.add_argument("-i", type=argparse.FileType("r"), default=sys.stdin,
      help="The input fasta file (default: stdin)")
parser.add_argument("-o", type=argparse.FileType("w"), default=sys.stdout,
      help="The output fasta file (default: stdout)")

args = parser.parse_args()

# args.mapping_file is a csv with 2 columnes 'old-name,new-name' 
# The <old-name> in input fasta will be replaced by <new-name>_<old_name>
mapping = {}
for old_name, new_name in csv.reader(args.mapping_file, dialect=csv.excel):
    mapping[old_name] = new_name

for line in args.i:
    if line.startswith(">"):
        name = line[1:].strip()
        if name in mapping:
            line = '>' + mapping[name] + '\n'
    args.o.write(line)

