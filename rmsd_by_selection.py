#!/Users/martinmccullagh/anaconda/bin/python
##!/Library/Frameworks/Python.framework/Versions/Current/bin/python
#USAGE : 
import scipy
import sys
import os
import numpy
import math
import MDAnalysis
from MDAnalysis.analysis.align import *

n_selections = 10

# defining which residues are each base type
res_select = []
res_select.append("bynum 10:20")
res_select.append("bynum 10:20")
res_select.append("bynum 10:20")
res_select.append("bynum 10:20")
res_select.append("bynum 10:20")
res_select.append("bynum 10:20")

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, traj_file, ref_file
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='topfile':
				top_file = value
			elif option.lower()=='trajfile':
				traj_file = value
			elif option.lower()=='reffile':
				ref_file = value
			else :
				print "Option:", option, " is not recognized"

# Main Program

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Trajectory file:", traj_file
print "Reference file:", ref_file

# initiate coordinate universe
coord = MDAnalysis.Universe(top_file, traj_file)
ref = MDAnalysis.Universe(top_file, ref_file)

# select nucleic acid backbone for alignment
coord_backbone = coord.selectAtoms("nucleicbackbone")
ref_backbone = ref.selectAtoms("nucleicbackbone")
# translate reference coordinates to origin
ref.translates(-ref_backbone.centerOfMass)
ref0 = ref_backbone.coordinates()

# make selections to compute RMSD for
for i in range(n_selections):
	sel[i] = coord.selectAtoms(res_select[i])
	ref_sel[i] = ref.selectAtoms(res_select[i])

# open output files
out = open("rmsd.dat",'w')

# Loop through trajectory
for ts in coord.trajectory:

	# align to reference
	coord.translate(-coord_backbone.centerOfMass())
	R, rmsd = rotation_matrix(center_backbone.coordinates(), ref0)
	coord.atoms.rotate(R)

	# write time step to output file
	out.write("%10d " % (ts.frame))

	# loop through selections and compute RMSD
	for i in range(n_selections):
		R, rmsd = rotation_matrix(sel[i].coordinates(),ref_sel[i].coordinates())
		out.write("%10.6f" % (rmsd))

	out.write("\n")

out.close

