#generates Ramachandran plots from Protein Data Bank (PDB) files

import sys, os
import numpy as np
import math
import matplotlib.pyplot as plt
import requests
import urllib

#Download the pdb through https GET request if not in current working directory
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

#Read pdb file, returns two arrays, groups to calculate phi and psi
def read_pdb(pdb_file):
	if not os.path.isfile(pdb_file + ".pdb"):
		print("Unable to find specified PDB file, downloading from RCSB...")
		download_pdb(pdb_file)

	bb_atoms = []
	phi_groups = []
	psi_groups = []
	bb_atom_names = ["N", "CA", "C"]
	chain = ""
	with open(pdb_file + ".pdb") as f:
		for line in f.read().splitlines():
			#Record backbone atom coordinates
			if line[0:4] == "ATOM" and line[13:17].strip() in bb_atom_names:
				#restarts at new chains - need?
				new_chain = line[20:23].strip()
				if new_chain != chain:
					bb_atoms = []
					chain = new_chain
				x = float(line[31:39])
				y = float(line[39:47])
				z = float(line[47:54])
				bb_atoms.append((x, y, z))

				atom_name = line[13:17].strip()
				#Group last four added and append to corresponding group
				if len(bb_atoms) > 3:
					if atom_name.strip() == "C":
						phi_groups.append(np.array(bb_atoms[-4:]))
						#phi_groups.append(bb_atoms[-4:])
					if atom_name.strip() == "N":
						psi_groups.append(np.array(bb_atoms[-4:]))
						#psi_groups.append(bb_atoms[-4:])

	return phi_groups, psi_groups

#Calculates dihedral given four coordinates (tuples of 3 floats)
def calculate_angle(a1, a2, a3, a4):
	#obtain the vectors b1, b2 and b3 by vector subtraction
	b1 = subtract(a2, a1)
	b2 = subtract(a3, a2)
	b3 = subtract(a4, a3)

	#Compute n1 = <b1 x b2> and n2 = ⟨<b2 x b3>
	#The angle we’re looking for is the same as the angle between n1 and n2
	n1 = normalize(cross_product(b1, b2))
	n2 = normalize(cross_product(b2, b3))

	#m1 = n1 x <b2>
	m1 = cross_product(n1, normalize(b2))

	#Compute the coordinates of n2 in this frame: x = n1 . n2 and y = m1 . n2
	x = dot_product(n1, n2)
	y = dot_product(m1, n2)

	#The dihedral angle, with the correct sign, is − atan2(y, x)
	return -math.atan2(y, x) * 180/math.pi

#Vector operations
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

#iterate through groups
def calculate_angles(coordinates):
	angles = []
	for coord in coordinates:
		angles.append(calculate_angle(*coord))
	return angles


#draw graph, could be prettier probably
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
#	3. Make a combined plot of 10,000 structures - done
#	4. Automatically download PDB files
#
##############################



