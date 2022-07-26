#!/usr/bin/env python
"""
Clean up ISOGG input from a .tsv

-Remove lines with non-canonical bases
-Remove lines with more than one mutation listed (ex: A->G; G->A)
-Remove sites that aren't single nucleotide changes (no indels)
-Remove sites with haplogroup names that end with "~", or that snp names that end with "^" 

For use with downstream script, fields are re-named and re-ordered according to 
Name    Subgroup_Name   Build_37_Number Build_38_Number Mutation_Info   Alternate_Names rs_number

This is an alpha release and may require spot-checking to ensure filtering is correct, 
and may require some additional curation

#need to check for if snps have more than one position for b37 and b38
#need to check for usable haplogroup names, possibly from reading in a tree?
## for example, one haplogroup shows as "R (notes)", which won't be helpful without cleaning
"""

__author__ = 'jashortt'

import argparse
import random
import os
import sys

def checkHeader(header):
	if (len(header) != 7):
		sys.exit("I was expecting 7 columns in the input but found %s." % (str(len(header))) )

def checkSNPs(snps):
	nucs = {'A', 'T', 'C', 'G'}
	keep_snps = 1
	if ( len(snps.split('->')) != 2 ):
		return(0)
	a1,a2= snps.split('->')
	for nuc in (a1,a2):
		if not nuc in nucs:
			return(0)
	return keep_snps
	
def checkHapName(hap):
	bad_chars = ["~", " "]
	keep_hap = 1
	if any(x in hap for x in bad_chars):	
		keep_hap = 0
	return keep_hap
	
def checkName(snp):
	keep_name = 1
	if("^" in snp):
		keep_name=0
	return keep_name
	
def clean_isogg (args):	
	outfile = "%s.txt" % (args.out)
	print ("Printing list of acceptable SNPS to %s." % (outfile))
	with open(args.infile, 'r') as infp, open(outfile, 'w') as outfp:
		garbage=infp.readline()
		header = infp.readline().strip().split('\t')
		checkHeader(header)
		#outfp.write( "%s\n"  %('\t'.join(x for x in header)) )
		outfp.write("Name\tSubgroup_Name\tBuild_37_Number\tBuild_38_Number\tMutation_Info\tAlternate_Names\trs_number\n")
		for line in infp:
			info = line.strip().split('\t')
			if (len(info) != 7):
				continue
			snp, hap, alt_snp, rs, pos37, pos38, snps = info
			keep_snps=checkSNPs(snps)
			keep_hap=checkHapName(hap)
			keep_name=checkName(snp)
			if (keep_snps == 1 & keep_hap == 1 & keep_name == 1):
				outfp.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snp, hap, pos37, pos38, snps, alt_snp, rs ))
			
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='cleanISOGG')

	parser.add_argument('--infile', help='ISOGG input, .tsv required', nargs='?', const=1, type=str, required=True)
	parser.add_argument('--out', help='prefix for output', nargs='?', const=1, type=str, default='isogg_snps', required=False)
	parser.add_argument('--version', action='version', version='%(prog)s alpha')
	
	args = parser.parse_args()
	clean_isogg(args)
	sys.exit()