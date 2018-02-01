#! /usr/bin/env python

import sys
import argparse
import requests,json
import re

def getOTTIDs(speciesString):
	'''
	Input is a string formatted with all the species names.
	Returns a response from the Open Tree of Life database with the matched species and their corresponding data and OTT ids.
	'''
	url = 'https://api.opentreeoflife.org/v3/tnrs/match_names' 
	headers = {'content-type':'application/json','Accept-Charset':'UTF-8'}
	data= '{"names":[' + speciesString + ']}' 
	response = requests.post(url,headers=headers,data=data).text
	return response

def readJSONresponse(response,matches,allOttIDs):
	'''
	Input is the JSON response from the OTL database.
	Also updates and returns  matches (a list of matched species from the response) and allOttIDs (a list of OTT ids found)
	'''
	data = json.loads(response)
	for result in data["results"]:
		ott_id = result["matches"][0]["taxon"]["ott_id"]
		numMatches = len(result["matches"])
		if numMatches == 1:
			allOttIDs.append(ott_id)
		elif result["matches"][0]["search_string"].lower() == result["matches"][0]["taxon"]["unique_name"].lower():
			allOttIDs.append(ott_id)
		else:
			matchList = []
			for match in result["matches"]:
				matchList.append(match)
			matches.append(matchList)
	return [matches,allOttIDs]

def removeDuplicates(matches,allOttIDs):
	'''
	Input is the list of matched species and list of all OTT ids
	Allows User to remove any duplicated species names
	Returns updated list of all OTT ids
	'''
	index = 0
	while index < len(matches):
		matchList = matches[index]
		search_string = matchList[0]["search_string"]
		ott_id = matchList[0]["taxon"]["ott_id"]
		print "Found multiple results for " + matchList[0]["matched_name"] + ". Please select the number for the species you want to keep."
		print 'Type "u" to undo.'
		count = 1
		matchesDict = {}
		for match in matchList:
			print str(count) + ". " + match["taxon"]["unique_name"]
			matchesDict[count] = match["taxon"]["unique_name"]
			count += 1
			
		userInput = getUserInput(index,count)
		if userInput == "u":
			allOttIDs = allOttIDs[:-1]
			index = index - 1
			continue
		val = int(userInput)
		ott_id = matchList[val - 1]["taxon"]["ott_id"]
		allOttIDs.append(ott_id)
		index += 1
	return allOttIDs

def getUserInput(index,count):
	'''
	Helps the removeDuplicates function get valid user input.
	Input is the current index number for the list of matches and the number of duplicate species found in the match.
	Returns the valid user input.
	'''
	badInput = True
	while badInput == True:
		answer = raw_input()
		if answer == "u":
			if index != 0:
				return "u"
			else:
				print "Sorry, nothing to undo. Please select a species."
		else:
			try:
				val = int(answer)
				if val < count and val > 0:
					badInput = False
				else:
					print "Invalid input. Please enter the number of your choice."
			except ValueError:
				print "Invalid input. Please enter the number of your choice."

	return answer

def makeNewickTree(allOttIDs):
	'''
	Input is our list of all OTT Ids
	Returns a response with the Newick tree from the OTL database.
	'''
	if len(allOttIDs) < 3:
		return False

	url ='https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree'
	headers = {'content-type':'application/json','Accept-Charset':'UTF-8'}
	data = '{"ott_ids":' + str(allOttIDs) + '}'
	response = requests.post(url,headers=headers,data=data).text

	if "The following OTT ids were not found: [" in response:  ####Remove id's that weren't found
		notFound = map(int,response.split("The following OTT ids were not found: [")[1].split("]. ")[0].encode('ascii','ignore').split(","))
		allOttIDs = list(set(allOttIDs) - set(notFound))
		if len(allOttIDs) < 3:
			return False
		data = '{"ott_ids":' + str(allOttIDs) + '}'
		response = requests.post(url,headers=headers,data=data).text
	newick = re.sub("_ott\d+","",response.split("newick\" : \"")[1].split(";\"")[0]) + ";"
	return newick

def printSpeciesNotFound(newick,specieslist,outfile):
	'''
	Input is the reference newick tree, the list of species from the original file, and the output file.
	Writes to the output file a list of all species that were not included in the reference tree.
	'''
	foundlist = re.split(',|\(|\)',newick.encode('ascii','ignore'))
	foundlist = [x for x in foundlist if len(x) > 0]
	foundlist = [x.replace("_"," ") for x in foundlist]
	allSpeciesNotFound = list(set(specieslist) - set(foundlist))
	if outfile != "":
		if len(allSpeciesNotFound) > 0:
			outfile.write("\nSpecies not found: " + ",".join(allSpeciesNotFound) + "\n")
	else:
		if len(allSpeciesNotFound) > 0:
			sys.stdout.write("\nSpecies not found: " + ",".join(allSpeciesNotFound) + "\n")

def readFile(args):
	'''
	Input is the system arguments.
	Reads through each line of the input file, makes a tree for each line, and outputs each line to its own output file.
	'''
	infile = open(args.input,"r")

	linenumber = 0
	specieslist = []
	for line in infile:
		linenumber += 1
		line = line.strip().replace("(","").replace(";","").replace(")","").replace("_"," ")
		specieslist.extend(line.split(","))
	if (len(specieslist) <= 2):
		print "Sorry, we didn't find enough species on line " + str(linenumber) + " to make a tree."
		return
	allOttIDs = formatOTTidRequests(specieslist)
	newick = makeNewickTree(allOttIDs)
	if newick == False:
		print "Sorry, we didn't find enough species on line " + str(linenumber) + " to make a tree."
		return
	newick = re.sub("\)[A-z_]+",")",newick)
	newick = re.sub("_"," ",newick)
	outfile = sys.stdout
	if args.output:
		outfile = open(args.output, "w")
	outfile.write(newick + "\n")
	if not args.excludeSpeciesNotFound:
		printSpeciesNotFound(newick,specieslist,outfile)
	if args.output:
		outfile.close()

def formatOTTidRequests(specieslist):
	'''
	Input is the list of species for which to make reference tree.
	Loops through the species 1000 at a time to retrieve OTT ids.
	Returns a list of all OTT ids.
	'''
	allOttIDs = []
	matches = []
	i = 0
	maxValue = 0
	while maxValue  < len(specieslist):
		minValue = i * 1000
		maxValue = (i + 1) * 1000
		if maxValue > len(specieslist):
			maxValue = len(specieslist)
		tempList = specieslist[minValue:maxValue]
		speciesString = ""
		if len(tempList) == 1:
			speciesString = tempList[0]
		else:	
			for species in tempList[:-1]:
				speciesString += "\"" + species + "\"" + ","
			speciesString += "\"" + tempList[-1] + "\""
		response = getOTTIDs(speciesString)
		i += 1
		if "\"matched_names\" : [ ]," in response:
			continue
		
		jsonResult = readJSONresponse(response,matches,allOttIDs)
		matches = jsonResult[0]
		allOttIDs = jsonResult[1]
	
	allOttIDs = removeDuplicates(matches,allOttIDs)
	return allOttIDs

def parseArgs():
	'''
	Does argument parsing.
	Required to have an input file.
	'''
	parser = argparse.ArgumentParser(description='Make Open Tree of Life Reference Tree from Newick Tree or Species List.')
	parser.add_argument("-i",help="Input Newick Tree or Species List",action="store", dest="input", required=True)
	parser.add_argument("-o",help="Output File Name",action="store",dest="output", required=False)
	parser.add_argument("-e",help="Exclude list of species not found in output file",action="store_true",dest="excludeSpeciesNotFound",required=False)
	args = parser.parse_args()
	
	return args


if __name__ == '__main__':
	'''
	Main.
	'''
	
	args = parseArgs()
	readFile(args)
