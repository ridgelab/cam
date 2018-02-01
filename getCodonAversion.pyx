#! /usr/bin/env python

import re


def getCodonAversion(sequence, possibleCodons):
	'''
	Takes two arguments: A DNA or RNA sequence, and a set of all possible codons.
	Subtracts a set of all codons in the sequence from the set of all possible codons.
	Returns a tuple of the subtracted set. This tuple represents all codons not found in the sequence.
	'''
	foundCodons= set(re.findall("...",sequence))
	motif = tuple(possibleCodons - foundCodons)
	return tuple(motif)
