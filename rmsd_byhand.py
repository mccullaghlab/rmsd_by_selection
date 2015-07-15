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

n_selections = 6

# defining which residues are each base type
res_select = []
res_select.append("bynum 3961:3965")
res_select.append("bynum 3994:3998")
res_select.append("bynum 4027:4031")
res_select.append("bynum 4059:4063")
res_select.append("bynum 4089:4093")
res_select.append("bynum 4121:4125")

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, traj_file, ref_file, out_file
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
			elif option.lower()=='outfile':
				out_file = value
			else :
				print "Option:", option, " is not recognized"

def computePbcDist2(r1,r2,box):
	dist2 = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist2 += temp*temp

	return dist2;

# Main Program

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Trajectory file:", traj_file
print "Reference file:", ref_file
print "Output RMSD data file:", out_file

# initiate coordinate universe
coord = MDAnalysis.Universe(top_file, traj_file)
ref = MDAnalysis.Universe(top_file, ref_file)

# select nucleic acid backbone for alignment
coord_backbone = coord.selectAtoms("nucleicbackbone")
ref_backbone = ref.selectAtoms("nucleicbackbone")
# translate reference coordinates to origin
ref.atoms.translate(-ref_backbone.centerOfMass())
ref0 = ref_backbone.coordinates()

# open output files
out = open(out_file,'w')

# Loop through trajectory
for ts in coord.trajectory:

	# align to reference
	coord.atoms.translate(-coord_backbone.centerOfMass())
	R, rmsd = rotation_matrix(coord_backbone.coordinates(), ref0)
	coord.atoms.rotate(R)

	# write time step to output file
	out.write("%10d " % (ts.frame))

	# loop through selections and compute RMSD
	for i in range(n_selections):
		sel = coord.selectAtoms(res_select[i])
		ref_sel = ref.selectAtoms(res_select[i])
		sel_rmsd = 0
		for atom in range(sel.numberOfAtoms()):
			sel_rmsd += computePbcDist2(sel.coordinates()[atom],ref_sel.coordinates()[atom],coord.dimensions[0:3])
		sel_rmsd = math.sqrt(sel_rmsd/float(sel.numberOfAtoms()))
		out.write("%10.6f" % (sel_rmsd))

	out.write("\n")

out.close

