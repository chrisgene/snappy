#!/usr/bin/env python
"""
create a list of SNPs to be used by SNAPPY
aims to have a user-specified count of snps for each haplotype
aims to cover every node from the root down to a user-specified distance (in node count) from the root

input is a .txt created from isogg_snp_compiler.py
"""

__author__ = 'jashortt'

import argparse
import sys

#returns pointer to a dictionary if key already exists or an empty dictionary
def dictIfEmpty(mydict, mykey):
	if not mykey in mydict:
		mydict[mykey]={}
	return(mydict[mykey])

#returns pointer to a set if key already exists or an empty set	
def setIfEmpty(mydict, mykey):
	if not mykey in mydict:
		mydict[mykey]=set()
	return(mydict[mykey])

#finds the name of the parent of a haplogroup	
def get_parent_hg(hg_to_parent, hg):
	"""returns the parental haplogroup"""
	if hg in hg_to_parent:
		return hg_to_parent[hg]
	elif hg in ['A0-T', 'A00']:
		return ''
	else:
		return hg[:-1]
  
#makes a list containing a haploroups ancestors back to the root          
def get_ancestry(hg, parent_dict):		
    """create list of all parent haplogroups"""
    ancestors = []
    ancestors.append(hg)
    hg = get_parent_hg(parent_dict, hg)
    while hg:
        ancestors.append(hg)
        hg = get_parent_hg(parent_dict, hg)
    return ancestors
    
#makes a dictionary of parent-child haplogroups that do not follow the standard nomenclature
def getHapTree (infile):
	group_to_parent = dict()
	with open(infile, 'r') as weird_hgs:
		for line in weird_hgs:
			line = line.rstrip('\n').split('\t')
			group_to_parent[line[0]] = line[1]
		return group_to_parent

#Creates a dictionary from input lines of a file. Each value is keyed by values in a header
def getLineInfo(line, header):
	info = dict()
	if len(header) != len(line):
		print('Incompatible line:\n%s\n%s' % ('\t'.join(header), '\t'.join(line)) )
		sys.exit()
	for key,val in zip(header, line):
		info[key] = val
	return info

#determines if a given mutation has strand-ambiguous alleles
def mutationIsAmbiguous (a1, a2):
	ambig_poly = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	if ambig_poly[a1] == a2:
		return 1
	else:
		return 0

#def determines if the alleles of a given mutation are transversions vs transitions. 
def mutationIsTransversion (a1, a2):
	good_poly = {'A':'G', 'C':'T', 'G':'A', 'T':'C'}			#might be more efficient to pass in than re-declare every time
	if good_poly[a1] != a2:
		return 1
	else:
		return 0
	
#determines mutation type for a snp		
def classifyMutationType(hap_info, info):
	snp_id = info['Name']
	if mutationIsAmbiguous(info['Ancestral_Allele'], info['Derived_Allele']):
		hap_info['ambig'] = dictIfEmpty(hap_info, 'ambig')
		hap_info['ambig'][snp_id] = info
	elif mutationIsTransversion(info['Ancestral_Allele'], info['Derived_Allele']):
		hap_info['trv'] = dictIfEmpty(hap_info, 'trv')
		hap_info['trv'][snp_id] = info
	else:						#mutation is a transition
		hap_info['trs'] = dictIfEmpty(hap_info, 'trs')
		hap_info['trs'][snp_id] = info									

#reads snp info from ISOGG file
def getHapSnpInfo(infile, mypos, build):
	print ('Getting SNP info from %s' % (infile))
	haps = {}
	hap_order = []
	with open(infile, 'r') as infp:
		header = infp.readline().strip().split()
		for line in infp:
			info = getLineInfo(line.strip().split('\t'), header)
			var_pos = info[build]
			if var_pos in mypos:				#only keep if variant is in position list
				hap = info['Subgroup_Name']				
				if not hap in haps:
					hap_info = dictIfEmpty(haps, hap)
					hap_order.append(hap)
				classifyMutationType(dictIfEmpty(haps, hap), info)
	return haps, hap_order
	
#makes a dictionary containing all ancestors for every haplogroup	
def getHapAncestry(haps, hap_order, hap_tree):
	for hap in hap_order:
		hap_info = dictIfEmpty(haps, hap)	
		hap_info['ancestors'] = get_ancestry(hap, hap_tree)

#gets build from args
def getHgBuild(build):
	pos_version = 'Build_37_Number'
	if build == 'hg38':
		pos_version ='Build_38_Number'
	return pos_version

#returns a set containing positions that were genotyped in a dataset
#the alleles SNAPPY chooses to make its reference library are based off these sites
def getPos(infile):
	mypos=set()
	with open(infile, 'r') as infp:
		for line in infp:
			var_pos = line.strip().split()
			if len(var_pos) > 1:
				print("Too many entries for this line: %s" % (" ".join(var_pos)))
				sys.exit()
			mypos.add(var_pos[0])
	return mypos
		
#select up to max_snp_count haplogroup informative snps for a haplogroup while prioritizing transition mutations over transversions, and transversions before strand-ambiguous snps	
def getSnps(hap_info, max_snp_count):
	hap_snps = []
	for snp_type in ['trs', 'trv', 'ambig']:
		if snp_type in hap_info:
			for snp in hap_info[snp_type].keys():
				snp_info = hap_info[snp_type][snp]
				hap_snps.append(snp_info)
	if len(hap_snps) <= max_snp_count:
		return hap_snps
	else:
		return hap_snps[0:max_snp_count-1]	

#get snps for each haplogroup
def makeSnps (haps, hap_order, max_snp_count, max_node_dist):
	snappy_snps = {}
	for hap in hap_order: #
		hap_info = haps[hap]
		if len(hap_info['ancestors']) <= int(max_node_dist): #
			snappy_snps[hap] = getSnps(hap_info, max_snp_count)
			keep_snps = snappy_snps[hap]
			#if len(keep_snps) < max_snp_count:
			#	print('Warning: only %s snps are available for %s' % (str(len(keep_snps)), hap))
			if len(keep_snps) == 0:
				print ('Warning: No snps found for %s' % (hap))
	return snappy_snps

#print snappy's reference files	
def printSnappyRefs (snappy_snps, hap_order, outfile, mybuild):	#might be nice to print a summary of snps and nodes in tree
	total_snp_count = 0
	hap_count = 0
	print ('Now opening %s to print snp list' % (outfile))
	with open(outfile, 'w') as outfp, open('id_to_pos.txt', 'w') as id2pos, open('pos_to_allele.txt', 'w') as pos2alleles, open('y_hg_and_snps.sort', 'w') as hgsort:
		id2pos.write('id\tpos\n')
		pos2alleles.write('pos\tancestral_allele\tderived_allele\n')
		hgsort.write('#haplogroup\tSNPs\n')
		for hap in hap_order:
			if hap in snappy_snps:
				hap_count += 1
				hap_snps = snappy_snps[hap]
				hgsort.write('%s\t%s\n' %( hap, ','.join([snp_info['Name'] for snp_info in hap_snps]) ))
				for snp_info in hap_snps:
					total_snp_count += 1
					snp_id = snp_info['Name']
					pos = snp_info[mybuild]
					ancestral = snp_info['Ancestral_Allele']
					derived = snp_info['Derived_Allele']
					outfp.write( '%s\n' % ( '\t'.join([hap, snp_id, pos, ancestral, derived])) )
					pos2alleles.write('%s\n' % ('\t'.join([pos, ancestral, derived])) )
					id2pos.write('%s\t%s\n' %(snp_id, pos))
	print ('Kept a total of total of %s snps for the %s haplogroups that met distance criterion.' % (str(total_snp_count), str(hap_count)))
				
def make_snappy_refs (args):
	mypos = getPos(args.pos_file)
	mybuild = getHgBuild(args.build)
	haps, hap_order = getHapSnpInfo(args.snp_file, mypos, mybuild)
	hap_tree = getHapTree(args.tree_file)
	getHapAncestry(haps, hap_order, hap_tree)
	keep_snps = makeSnps(haps, hap_order, args.snp_count, args.max_node_dist)
	printSnappyRefs	(keep_snps, hap_order, '%s.txt' % (args.out), mybuild)
		
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='make_snappy_refs', description="Make a set of reference files for use by SNAPPY. Makes some attempt at providing balanced allelic representation across Y-chromosome phylogeny but this is inherently limited by the SNPs present in the input data")
    
    parser.add_argument('--snp_file', help='a tab-separated file containing ISOGG snps', required=True)
    parser.add_argument('--pos_file', help='list of physical positions of Y-chromosome genotypes, no header', required=True)
    parser.add_argument('--tree_file', help='list of non-canonical relationships between haplogroups names, distributed with SNAPPY', default='ref_files/tree_structure.txt', required=False)
    parser.add_argument('--snp_count', help='a target for the number of haplogroup-informative SNPs to include for each haplogroup', type=int, default=5, required=False)
    parser.add_argument('--max_node_dist', help='the maximim distance from the root to build the tree', type=int, default=99, required=False)
    parser.add_argument('--build', help='genome build, hg37 or hg38', nargs='?', const='hg37', choices=['hg37', 'hg38'], type=str, default='hg37', required=False)
    parser.add_argument('--out', help='prefix for file output', nargs='?', const=1, type=str, default='SNAPPY_snp_list', required=False)

    args = parser.parse_args()
    make_snappy_refs(args)
    sys.exit()