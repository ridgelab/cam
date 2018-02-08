#! /usr/bin/env python
'''
Reads a distance matrix file created by cam.py and
writes a newick tree to standard out.

'''
import sys
import os
import argparse

def checkTempNum(TEMP_FILE_NUM):
	'''
	Ensures that the same temporary file is not used.
	'''
	for fname in os.listdir('.'):
	    if fname.endswith(TEMP_FILE_NUM):
			return True
	return False


def getSpeciesDistances(fileName):
	'''
	Input is the path to the distance matrix.
	Returns a list of species names and the top half of the distance matrix.
	'''
	inputFile = open(fileName,'r')
	species = inputFile.readline().strip().split(',')[1:]
	distance = []
	pos = 2 #Start at position 2 because position 0 is the species name and position 1 is always 0 (distance to itself).
	for line in inputFile:
		distance.append(map(float,line.strip().split(",")[1:pos]))
		pos +=1 ####Commented out when only top of matrix present
	return species, distance

def writeNewick(species, distance,output):
	'''
	Input is a list of species names and the top half of the distance matrix
	Newick tree is output to standard out.	

	'''
	outputFile = sys.stdout
	if args.output:
		outputFile = open(output,'w')
	import Bio.Phylo.TreeConstruction  as TreeConstruction
	constructor = TreeConstruction.DistanceTreeConstructor()
	distanceMatrix = TreeConstruction._DistanceMatrix(species,distance)
	treeConstructor = TreeConstruction.DistanceTreeConstructor(method = 'nj')
	njTree = treeConstructor.nj(distanceMatrix)
	TEMP_FILE_NUM = str(int(os.urandom(3).encode('hex'),16))
	while checkTempNum(TEMP_FILE_NUM):
		TEMP_FILE_NUM = str(int(os.urandom(3).encode('hex'),16))
	tempFile = open(".tempFile" + TEMP_FILE_NUM,'w')	
	from Bio import Phylo
	Phylo.write(njTree,tempFile,"newick")
	tempFile.close()
	import re
	treeF = open(".tempFile" + TEMP_FILE_NUM,'r')	
	tree = treeF.read()
	treeF.close()
	os.remove(".tempFile" +TEMP_FILE_NUM)

	tree = re.sub("Inner[0-9]+:[-0-9\.]+","",tree)
	tree = re.sub(":[0-9\.]+","",tree)
	tree = re.sub("_"," ",tree)
	outputFile.write(tree)
	if args.output:
		outputFile.close()

def parseArgs():
	'''
	Argument parsing is done.
	Required to have an input file.
	'''
	parser = argparse.ArgumentParser(description='Make Newick File from Distance Matrix.')
	parser.add_argument("-i",help="Input Fasta Files",action="store", dest="input", required=True)
	parser.add_argument("-o",help="Output File",action="store",dest="output", required=False)
	args = parser.parse_args()
	
	return args



if __name__ =='__main__':
	'''
	Main.
	'''
	args = parseArgs()
	fileName = sys.argv[1]
	species,distance = getSpeciesDistances(args.input)
	writeNewick(species,distance,args.output)
