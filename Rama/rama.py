#generates Ramachandran plots from Protein Data Bank (PDB) files

import sys, os
import numpy as np
import math

def read_pdb(pdb_file):
	if not os.path.isfile(pdb_file):
		print("ERROR: Unable to find specified PDB file")
		sys.exit(-1)

	bb_atoms = []
	phi_groups = []
	psi_groups = []
	bb_atom_names = ["CA", "C", "N"]
	with open(pdb_file) as f:
		for line in f.read().splitlines():
			if line[0:4] == "ATOM" and line[13:17].strip() in bb_atom_names:
				resn = int(line[23:27])
				x = float(line[31:39])
				y = float(line[39:47])
				z = float(line[47:54])
				bb_atoms.append((x, y, z))

				if len(bb_atoms) > 3:
					if line[13:17].strip() == "C":
						phi_groups.append(np.array(bb_atoms[-4:]))
					if line[13:17].strip() == "N":
						psi_groups.append(np.array(bb_atoms[-4:]))

	return phi_groups, psi_groups

def calculate_angle(a1, a2, a3, a4):
	b1 = subtract(a1, a2)
	b2 = subtract(a2, a3)
	b3 = subtract(a3, a4)

	n1 = cross_product(b1, b2)/normalize(cross_product(b1, b2))
	n2 = cross_product(b2, b3)/normalize(cross_product(b2, b3))

	m1 = cross_product(n1, b2/normalize(b2))

	x = dot_product(n1, n2)
	y = dot_product(m1, n2)

	return math.atan2(y, x) * 180/math.pi


def add(v1, v2):
	return np.add(v1, v2)

def subtract(v1, v2):
	return np.subtract(v1, v2)

def normalize(v1):
	return np.linalg.norm(v1)

def dot_product(v1, v2):
	return np.dot(v1, v2)

def cross_product(v1, v2):
	return np.cross(v1, v2)

def calculate_angles(coordinates):
	angles = []
	for coord in coordinates:
		angles.append(calculate_angle(*coord))
	return angles

def main():
	if len(sys.argv[1:]) < 1:
		print("ERROR: No PDB file specified")
		sys.exit(-1)
	phis, psis = read_pdb(sys.argv[1])
	phi_angles = calculate_angles(phis)
	print(phi_angles)


if __name__ == "__main__":
	main()