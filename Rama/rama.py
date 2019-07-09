#generates Ramachandran plots from Protein Data Bank (PDB) files

import sys


def main():
	if len(sys.argv[1:]) < 1:
		print("Error: No PDB file specified")
		sys.exit()
	make_ramachandran(sys.argv[1])