#! /usr/bin/env python
import sys
import argparse
from getCodonAversion import getCodonAversion
from multiprocessing import Process, current_process, freeze_support, Pool
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna


def makeAllPossibleCodons(args):
	'''
	Input is a flag to specify if the sequence is DNA or RNA.
	Returns a set of all 64 possible codons.
	'''
	if args.amino or args.aa:
		return set(['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
	from itertools import product
	codons = product("ACGT",repeat=3)
	if args.rna:
		codons = product("ACGU",repeat=3)
	codonsComb = set()
	for c in codons:
		codonsComb.add("".join(c))
	return codonsComb

def parseArgs():
	'''
	Argument parsing is done.
	Required to have an input file.
	'''
	parser = argparse.ArgumentParser(description='Make Distance Matrix from Codon Aversion Motifs.')
	parser.add_argument("-t",help="Number of Cores",action="store",dest="threads",default=0,type=int, required=False)
	parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
	parser.add_argument("-id",help="Input Directory with Fasta Files",action="store", dest="inputDir", required=False)
	parser.add_argument("-o",help="Output File",action="store",dest="output", required=False)
	parser.add_argument("-p",help="Percent of codon aversion tuples Shared By Species",action="store",dest="percent", type=int, default=0.05, required=False)
	parser.add_argument("-rna",help="Flag for RNA sequences",action="store_true",dest="rna", required=False)
	parser.add_argument("-a",help="Flag for amino acid sequences",action="store_true",dest="amino", required=False)
	parser.add_argument("-aa",help="Flag for Amino Acid Motifs (From DNA/RNA)",action="store_true",dest="aa", required=False)
	parser.add_argument("-w",help="Writes the distance matrix immediately to the output file. Header will be at the bottom.",action="store_false",dest="write", required=False) 
	args = parser.parse_args()
	
	if not args.input and not args.inputDir:
		print "You must supply an input file with either -i or -id"
		sys.exit()
	return args

def readOneFile(inputFile):
	'''
	Reads one input file that is supplied as a parameter.
	Returns a tuple of the input file name and the set of tuples where each
		tuple is a set of codons not present in at least one gene.
	'''
	input = ""
	setOfTuples = set()
	header = ""
	sequence = ""
	try:
		if inputFile[-3:] =='.gz':
			import gzip
			input = gzip.open(inputFile,'r')
		else:
			input = open(inputFile,'r')
		for line in input:
			if line[0] =='>':
				if sequence !="":
					if args.aa:
						aa =""
						if args.rna:
							rna = Seq(sequence,generic_rna)
							aa = str(rna.translate())
						else:
							seq = Seq(sequence,generic_dna)
							aa = str(seq.translate())
						setOfTuples.add(tuple(set(list(aa))))
					elif args.amino:
						setOfTuples.add(tuple(set(list(sequence))))
					else:
						setOfTuples.add(getCodonAversion(sequence,codonsComb))
				header = line.strip()
				sequence = ""
				continue
			sequence +=line.upper().strip()
		if sequence != "":
			if args.aa:
				aa =""
				if args.rna:
					rna = Seq(sequence,generic_rna)
					aa = str(rna.translate())
				else:
					seq = Seq(sequence,generic_dna)
					aa = str(seq.translate())
				setOfTuples.add(tuple(set(list(aa))))
			elif args.amino:
				setOfTuples.add(tuple(set(list(sequence))))
			else:
				setOfTuples.add(getCodonAversion(sequence, codonsComb))
			sequence = ""
	except Exception: #If the input file is malformatted, do not stop the program.
		input.close()
		return set()
	input.close()
	return ((inputFile,setOfTuples))

def readInputFiles(args):
	'''
	Requires arguments to be passed to the function.
	Returns a dictionary where the keys are the input file names 
		and the values are sets of tuples where the tuples are codons not present in 
		at least one gene in the species file.
	'''
	threads = args.threads
	if threads ==0:
		pool = Pool()
	else:
		pool = Pool(threads)	
	allInputFiles = []
	allSets = set()
	fileToSet = {}
	if args.input:
		allInputFiles = args.input
	elif args.inputDir:
		import os
		allFasta = []
		path = args.inputDir
		allFasta = os.listdir(path)
		if path[-1] != '/':
			path += '/'
		allInputFiles = [path +i for i in allFasta]
	if len(allInputFiles) < 2:
		print "At least two input files are required"
		sys.exit()
	temp =pool.map(readOneFile,allInputFiles)
	for x in temp:
		if len(x) ==0:
			continue
		inputFile = x[0]
		setOfTuples = x[1]
		allSets |=setOfTuples
		fileToSet[inputFile] = setOfTuples
	if args.output:
		print "Total Number of motifs:",len(allSets)
	return fileToSet

def writeDistanceMatrix(args):
	'''
	Takes as input the system arguments
	Writes the Distance matrix to the output file or to standard out
	'''
	output = ""
	if args.output:
		output = open(args.output,'w')
	fileNames = []
	fullOutputString = ""
	for file1 in fileToSet:
		fileNames.append(file1.split("/")[-1].split(".gz")[0]) #Adds the name of the species to a list
		outString = fileNames[-1]
		for file2 in fileToSet:
			overlapping = fileToSet[file1] & fileToSet[file2]
			possible = float(min(len(fileToSet[file1]),len(fileToSet[file2])))
			maxPossible = float(max(len(fileToSet[file1]),len(fileToSet[file2])))
			if maxPossible == 0:
				outString += ",1"
				continue
			possiblePercent = possible/ maxPossible
			distance = possible
			included = float(len(overlapping))
			if possiblePercent > args.percent: 
				distance = 1 - (included /possible)
			else:
				distance = "1"
			outString += ("," + str(distance))
		outString += "\n"
		if args.write:
			fullOutputString += outString #Adds a row to the distance matrix
		else:
			if args.output:
				output.write(outString)
			else:
				sys.stdout.write(outString)
	s = "," + ",".join(fileNames) + "\n" #Creates the header line
	if args.write:
		if args.output:
			output.write(s + fullOutputString)
		else:
			sys.stdout.write(s + fullOutputString)
	else:
		if args.output:
			output.write(s)
			output.close()
		else:
			sys.stdout.write(s)

if __name__ =='__main__':
	'''
	Main.
	'''

	args = parseArgs()
	codonsComb = makeAllPossibleCodons(args)
	fileToSet = readInputFiles(args)
	writeDistanceMatrix(args)


