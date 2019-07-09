#generates Ramachandran plots from Protein Data Bank (PDB) files

import sys, os
import numpy as np
import math
import matplotlib.pyplot as plt
import requests
import urllib

def download_pdb(pdb_name):
	url = 'https://files.rcsb.org/download/' + pdb_name + '.pdb'
	r = requests.get(url)
	print(r.status_code)
	if r.status_code == 200:
		with open(pdb_name + '.pdb', "w+") as f:
			f.write(r.text)
	else:
		print("Error retrieving pdb file, likely due to incorrect PDB code.")
		sys.exit(-1)

def read_pdb(pdb_file):
	if not os.path.isfile(pdb_file + ".pdb"):
		print("Unable to find specified PDB file, downloading from RCSB...")
		download_pdb(pdb_file)

	bb_atoms = []
	phi_groups = []
	psi_groups = []
	bb_atom_names = ["CA", "C", "N"]
	chain = ""
	with open(pdb_file + ".pdb") as f:
		for line in f.read().splitlines():
			if line[0:4] == "ATOM" and line[13:17].strip() in bb_atom_names:
				if line[20:23].strip() != chain:
					bb_atoms = []
					chain = line[20:23].strip()
				resn = int(line[23:27])
				x = float(line[31:39])
				y = float(line[39:47])
				z = float(line[47:54])
				bb_atoms.append((x, y, z))

				if len(bb_atoms) > 3:
					if line[13:17].strip() == "C":
						phi_groups.append(np.array(bb_atoms[-4:]))
						#phi_groups.append(bb_atoms[-4:])
					if line[13:17].strip() == "N":
						psi_groups.append(np.array(bb_atoms[-4:]))
						#psi_groups.append(bb_atoms[-4:])

	return phi_groups, psi_groups


def calculate_angle(a1, a2, a3, a4):
	b1 = subtract(a1, a2)
	b2 = subtract(a2, a3)
	b3 = subtract(a3, a4)

	n1 = normalize(cross_product(b1, b2))
	n2 = normalize(cross_product(b2, b3))

	m1 = cross_product(n1, normalize(b2))

	x = dot_product(n1, n2)
	y = dot_product(m1, n2)

	return math.atan2(y, x) * 180/math.pi


def add(v1, v2):
	#return tuple(map(lambda x: sum(x), zip(v1, v2)))
	return np.add(v1, v2)

def subtract(v1, v2):
	#return tuple(map(lambda x: x[0] - x[1], zip(v1, v2)))
	return np.subtract(v1, v2)

def normalize(v1):
	#length = math.sqrt(sum(map(lambda x: x**2, v1)))
	#return tuple(map(lambda x: x/length, v1))
	length = np.linalg.norm(v1)
	return np.true_divide(v1, length)

def dot_product(v1, v2):
	#return sum(map(lambda x: x[0]*x[1], zip(v1, v2)))
	return np.dot(v1, v2)

#couldn't think of a fun way to do this one :(
def cross_product(v1, v2):
	#x = v1[1]*v2[2] - v1[2]*v2[1]
	#y = v1[2]*v2[0] - v1[0]*v2[2]
	#z = v1[0]*v2[1] - v1[1]*v2[0]

	#return (x, y, z)
	return np.cross(v1, v2)

def calculate_angles(coordinates):
	angles = []
	for coord in coordinates:
		angles.append(calculate_angle(*coord))
	return angles

def graph_ramachandran(name, phi_angles, psi_angles):
	fig, ax = plt.subplots()
	ax.scatter(phi_angles, psi_angles, marker='.')
	ax.set_xlim(-180, 180)
	ax.set_ylim(-180, 180)
	ax.set_adjustable('box')
	ax.set_aspect(1./ax.get_data_ratio())
	ax.set_xlabel("Φ")
	ax.set_ylabel("Ψ")
	fig.tight_layout()
	fig.savefig(name + '_rama.png', dpi=300)

if len(sys.argv[1:]) < 1:
	print("ERROR: No PDB code specified")
	sys.exit(-1)

name = sys.argv[1].replace(".pdb", "")
phis, psis = read_pdb(name)
phi_angles = calculate_angles(phis)
psi_angles = calculate_angles(psis)
graph_ramachandran(name, phi_angles, psi_angles)

### Programming challenges ###
#
#	1. Plot interesting subsets of residues
#	2. Use NumPy to speed up vector arithmetic - done
#	3. Make a combined plot of 10,000 structures
#	4. Automatically download PDB files
#
##############################



