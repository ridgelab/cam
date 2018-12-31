#! /usr/bin/env python
'''
Reads a distance matrix file created by cam.py and
writes a newick tree to standard out.

'''
import sys
import os
import argparse
from base64 import b64encode

def checkTempNum(TEMP_FILE_NUM):
    '''
    Ensures that the same temporary file is not used.
    '''
    for fname in os.listdir('.'):
        if fname.endswith(TEMP_FILE_NUM):
            return True
    return False

def getMin(distance,species):
    values = dict()
    for x in distance: ##x = species
        myMin = min(distance[x])
        minY = distance[x].index(myMin)
        minX = x
        if not myMin in values:
            values[myMin] = []
        values[myMin].append(tuple((minX,species[minY])))
    return values

def recalibrateDistance(distance,species,values):
    for op in sorted(values.keys()):
        if op==1:
            return distance,species
        for listZ in values[op]: 
            x,y = listZ
            if not x in species or not y in species:
                return distance,species
            indexX = species.index(x)
            indexY = species.index(y)
            speciesX = '(' + y +"," + x + ")"
            species[indexX] = speciesX
            distance[speciesX] = []
            for z in range(len(distance[x])):
                distance[speciesX].append(min(distance[x][z],distance[y][z]))
            del distance[x]
            del distance[y]
            for s in distance:
                del distance[s][indexY]

            del species[indexY]
    return distance,species



def getTree(species, distance): 
    lastLength = len(species)    
    while len(species)>1:
        sys.stderr.write("Species Remaining=" + str(len(species)) +"\n" )
        values =getMin(distance,species)
        distance,species = recalibrateDistance(distance,species,values)
        if len(species) == lastLength:
            break
        lastLength = len(species)
    return species

def largeSpeciesTree(args):
    fileNames = []
    matrix = dict()
    inputFile = open(args.input)
    output = sys.stdout
    if args.output:
        output = open(args.output,'w')
    fileNames = inputFile.readline().strip().split(",")[1:]
    numToDel = 0
    for line in inputFile:
        sys.stderr.write("Species Added to Matrix=" + str(len(matrix)) +"\n" )
        numToDel +=1
        distances = line.strip().split(",")
        species = distances[0]
        matrix[species] = []
        distances = map(float, distances[1:])
        num = 0
        for x in distances:
            if num<numToDel:
                matrix[species].append(1) #Max value is 1 in matrix
                num +=1
            else:
                matrix[species].append(x)
                continue
    distTree = getTree(fileNames,matrix)    
    newTree = ""
    if len(distTree)>1:
        newTree += '('
    for x in range(len(distTree)):
        newTree += distTree[x].replace('_',' ')
        if x < len(distTree)-1:
            newTree += ','
    if len(distTree)>1:
        newTree += ')'
    newTree += ';'
    output.write(newTree)
    output.close()

def getSpeciesDistances(args):
    '''
    Input is the path to the distance matrix.
    Returns a list of species names and the bottom half of the distance matrix.
    '''
    if args.phylip:
        species = []
        distance = []
        inputFile = open(args.input,'r')
        inputFile.readline()
        pos = 1 
        for line in inputFile:
            distance.append(list(map(float,line[10:].strip().split(" "))[0:pos]))
            species.append(line[0:10].strip())
            pos +=1
        return species,distance
            

    inputFile = open(args.input,'r')
    species = inputFile.readline().strip().split(',')[1:]
    if len(species) >500: #if it's a large tree
        inputFile.close()
        largeSpeciesTree(args)
        sys.exit()
    distance = []
    pos = 2 #Start at position 2 because position 0 is the species name and distance to itself (0.0) needs to be included.
    for line in inputFile:
        distance.append(list(map(float,line.strip().split(",")[1:pos])))
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
    TEMP_FILE_NUM = b64encode(os.urandom(3)).decode('utf-8')
    while checkTempNum(TEMP_FILE_NUM):
        TEMP_FILE_NUM = b64encode(os.urandom(3)).decode('utf-8')
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
    tree = re.sub(":[-0-9\.]+","",tree)
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
    parser.add_argument("-p",help="Phylip format",action="store_true",dest="phylip", required=False)
    args = parser.parse_args()
    
    return args



if __name__ =='__main__':
    '''
    Main.
    '''
    args = parseArgs()
    fileName = sys.argv[1]
    species,distance = getSpeciesDistances(args)
    writeNewick(species,distance,args.output)
