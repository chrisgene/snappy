#!/usr/bin/env python
"""
compile lists of snps for Y-chromosome haplogroups
input is a curated list from ISOGG (https://isogg.org/tree/index.html)
"""

__author__ = 'jashortt'

import argparse
import os
import sys
import re
from copy import deepcopy
from snappy import *

def dictIfEmpty(mydict, mykey):
	if not mykey in mydict:
		mydict[mykey]={}
	return(mydict[mykey])
	
def setIfEmpty(mydict, mykey):
	if not mykey in mydict:
		mydict[mykey]=set()
	return(mydict[mykey])
				
def getLineInfo(line, header):
	info = dict()
	for key,val in zip(header, line):
		info[key] = val
	return info
	
def getAlleles(mutation):
	''' returns list of alleles in following order: ancestral, derived'''
	myalleles = re.split('->|-|>', mutation)			#check that this is working correctly
	if len(myalleles) != 2:
		print('There can be only two alleles per line. Offending mutation info: in: %s\tout: %s' % (mutation, ' '.join(myalleles)))
		print ('Please clean input file before compiling snps.')
		sys.exit()
	for allele in myalleles:
		if not allele in ('A', 'T', 'G', 'C'):
			print ('Abnormal mutation info: in: %s\tout: %s\toffending allele: %s' % (mutation, ' '.join(myalleles), allele))
			print ('Please clean input file before compiling snps.')
			sys.exit()
	return {'ancestral': myalleles[0], 'derived': myalleles[1]}
	
def getAltNames(info):						#need to check that this split is working well
	alts = {info['Name']}
	if 'Alternate_Names' in info and info['Alternate_Names'] != '':
		alts.update(re.split(';|,', info['Alternate_Names']))
	return alts
					
def getHgBuild(build):
	pos_version = 'Build_37_Number'
	if build == 'hg38':
		pos_version ='Build_38_Number'
	return pos_version

def saveAlleleInfo(mydict, info):
	info['alts'] = getAltNames(info)
	var_id = info['Name']
	if not var_id in mydict:      #will only save info the first time this id shows up. If it shows up a second time with different info, the different info won't be saved
		mydict[var_id] = info
	for alt in info['alts']:
		if not alt in mydict:
			mydict[alt] = mydict[var_id]
		
def getInfoByPos(infile, mybuild):
	pos = {}
	with open(infile, 'r') as infp:
		header = infp.readline().strip().split()
		for line in infp:
			info = getLineInfo(line.strip().split('\t'), header)
			snp_pos = info[mybuild]
			info['alleles'] = getAlleles(info['Mutation_Info'])
			alt_names = getAltNames(info)
			pos_dict = dictIfEmpty(pos, snp_pos)
			anc_dict = dictIfEmpty(pos_dict, info['alleles']['ancestral'])
			der_dict = dictIfEmpty(anc_dict, info['alleles']['derived'])
			saveAlleleInfo(der_dict, info)
	return pos
	
def removeBadInfo(pos, derived_to_remove, ancestral_to_remove, pos_to_remove):
	for mypos in derived_to_remove.keys():
		for anc in derived_to_remove[mypos].keys():
			for der in derived_to_remove[mypos][anc]:
				del pos[mypos][anc][der]
	for mypos in ancestral_to_remove.keys():
		for anc in ancestral_to_remove[mypos].keys():
				del pos[mypos][anc]
	for mypos in pos_to_remove:
		del pos[mypos]
		
def haplogroupCheck(pos):
	'''
	Runs through each set of ancestral and derived alleles at each position and removes sets with more than one haplogroup listed for it
	'''
	cpos = deepcopy(pos)
	for snp_pos in cpos.keys():
		for anc in cpos[snp_pos].keys():
			for der in cpos[snp_pos][anc]:
				haps = set()
				for var_id in cpos[snp_pos][anc][der]:
					info = cpos[snp_pos][anc][der][var_id]
					haps.add(info['Subgroup_Name'])
				hap_count = len(haps)
				if(hap_count) > 1:
					del pos[snp_pos][anc][der]
			derived_count = len(pos[snp_pos][anc].keys())
			if derived_count < 1:
				del pos[snp_pos][anc]
		anc_count = len(pos[snp_pos].keys())
		if anc_count < 1:
			del pos[snp_pos]
	
def var_IdCheck(pos):
	var_ids = {}
	for snp_pos in pos.keys():
		for anc in pos[snp_pos].keys():
			for der in pos[snp_pos][anc]:
				info = pos[snp_pos][anc][der]
				for var_id in info.keys():
					id_info = dictIfEmpty(info, var_id)
					id_pos = setIfEmpty(id_info, 'pos')
					id_pos.add(snp_pos)
					id_haps = setIfEmpty(id_info, 'haps')
					id_haps.add(id_info['Subgroup_Name'])
	#check for multiple haps first, if multiple, remove id using 'pos'
	for var_id in var_ids.keys():
		hap_count = len(var_ids[var_id]['haps'])
		if hap_count > 1:
			for var_pos in var_ids[var_id]['pos']:
				del pos[var_pos]			#might be a little harsh to delete entire pos... other vars at site could be well-behaved
				del var_ids[var_id]
	for var_id in var_ids.keys():
		pos_count = len(var_ids[var_id]['pos'])
		if pos_count > 1:
			del var_ids[var_id]
			for var_pos in var_ids[var_id]['pos']:
				del pos[var_pos]		#again, might be too harsh to remove the entire pos, if it's just one id at the pos acting up			
	
def checkInfo(pos):
	'''
	Do some quality control on the data:
		a given derived allele at a position cannot have more than one haplogroup
		a given id should not appear at multiple positions or with multiple haplogroups
	'''
	pos_count = len(pos.keys())
	if pos_count < 1:
		print("There are no positions remaining before quality control. Goodbye")
		sys.exit()
	else:
		print("There are %s positions remaining before quality control." % (str(pos_count)))
	haplogroupCheck(pos)
	pos_count = len(pos.keys())
	if pos_count < 1:
		print("There are no positions remaining after haplogroup quality control. Goodbye")
		sys.exit()
	else:
		print("There are %s positions remaining after haplogroup quality control." % (str(pos_count)))
	#run through each position, save each varid and its position. Remove those that have more than one position and/or haplogroup
	var_IdCheck(pos)
	pos_count = len(pos.keys())
	if pos_count < 1:
		print("There are no positions remaining after variant id quality control. Goodbye")
		sys.exit()
	else:
		print("There are %s positions remaining after variant id quality control." % (str(pos_count)))
	
def printVars (outfile, pos):
	print("Printing list of quality SNPs to %s" %(outfile))
	with open(outfile, 'w') as outfp:
		outfp.write('Name\tSubgroup_Name\tAlternate_Names\tBuild_37_Number\tBuild_38_Number\tAncestral_Allele\tDerived_Allele\n')
		for snp_pos in sorted(pos.keys()):
			for anc in pos[snp_pos].keys():
				for der in pos[snp_pos][anc]:
					printed = set()
					allele_info = pos[snp_pos][anc][der]
					for var_id in allele_info:
						snp_info = pos[snp_pos][anc][der][var_id]
						if not snp_info['Name'] in printed:  #checking if an alternate name of the snp is already printed
							printed.add(var_id)
							snp_id = snp_info['Name']
							hap = snp_info['Subgroup_Name']
							alts = ""
							if 'alts' in snp_info:
								alts = ','.join( snp_info['alts'] )
							b37_pos = snp_info['Build_37_Number']
							b38_pos = snp_info['Build_38_Number']
							mutation = '\t'.join([snp_info['alleles']['ancestral'], snp_info['alleles']['derived']])
							outfp.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (snp_id, hap, alts, b37_pos, b38_pos, mutation) )
	
#make a haplogroup dictionary where info on haplogroup-informative alleles is easily accessed
def isogg_qc (args):
	mybuild = getHgBuild(args.build)
	pos = getInfoByPos(args.infile, mybuild)
	checkInfo(pos)
	printVars('%s.txt' % (args.out), pos)
		
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--infile', help='a tab-separated file containing ISOGG snps', required=True)
    parser.add_argument('--build', help='genome build, hg37 or hg38', nargs='?', const='hg37', choices=['hg37', 'hg38'], type=str, default='hg37', required=False)
    parser.add_argument('--out', help='prefix for file output', nargs='?', const='snp_qc', type=str, default='snp_qc', required=False)

    args = parser.parse_args()
    isogg_qc(args)
    sys.exit()
