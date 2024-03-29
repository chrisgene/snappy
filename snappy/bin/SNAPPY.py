"""
SNAPPY.py, by Alissa Severson, with Jonathan A. Shortt

This script takes chrY data in plink format (.bim, .bed, .fam) and produces haplogroup assignments.
"""

import os
import os.path
import sys
import subprocess
from snappy import *
from snappy.bin.parse_ref_files import *        #these two lines aren't clean-looking but they do at least seem to work
from snappy.bin.parse_plink_files import *
import argparse

def count_define_called_snps(subgroup_snps, pos_to_gt, isogg_id_to_pos):
    """counts how many defining snps are actually called"""
    total = 0
    for snp in subgroup_snps:
        # if multiple ids for the same snp, tests until finds one with position
        if '/' in snp:
            candidates = snp.split('/')
            for candidate in candidates:
                if candidate in isogg_id_to_pos:
                    snp = candidate
                    break
        # if the snp is genotyped, count it
        if snp in isogg_id_to_pos:
            pos = isogg_id_to_pos[snp]
            if pos in pos_to_gt:
                gt = pos_to_gt[pos]
                if gt != 'N':
                    total += 1
    return float(total)


def get_parent_hg(hg_to_parent, hg):
    """returns the parental haplogroup"""
    if hg in hg_to_parent:
        return hg_to_parent[hg]
    elif hg in ['A0-T', 'A00']:
        return ''
    else:
        return hg[:-1]


def has_parent_calls(hg, sample, scores, names, parent_dict, ancestral_hg_depth):
    """check if at least the parent or grandparent hg has a derived genotype"""
    count_present = 0   # keeps track of whether parent or grandparent hg is derived
    for i in range(ancestral_hg_depth):
        hg = get_parent_hg(parent_dict, hg)     # get parental hg
        if hg in names:
            hg_index = names.index(hg)
            if scores[sample, hg_index] > 0:    # check if hg snps are derived
                count_present += 1

    if count_present > 0:
        return True
    else:
        return False


def get_ancestry(hg, parent_dict):
    """create list of all parent haplogroups"""
    hg = get_parent_hg(parent_dict, hg)
    ancestors = []
    while hg:
        ancestors.append(hg)
        hg = get_parent_hg(parent_dict, hg)

    return ancestors


def score_hgs(hg_scores, hg_to_snps, genotypes, n, issog_id_to_pos, group_to_parent, ancestral_hg_depth):
    """
    score every hg for an individual using the counts recorded in hg_score, calculate what fraction of snps and
    ancestral snps are derived
    """
    hg_names = list(hg_to_snps.keys())
    strict_hg_to_score = dict()     # only score hgs with derived parent or grandparent hg
    all_hg_to_score = dict()        # score all hgs
    for h in range(len(hg_names)):
        hg_score = hg_scores[n, h]  # get number of derived calls for hg snps
        if hg_score > 0:
            hg = hg_names[h]
            n_called_snps = 0
            n_defining_snps = 0
            while hg:               # count number of derived snps in hg and all ancestral hgs
                if hg in hg_names:
                    n_called_snps += hg_scores[n, hg_names.index(hg)]
                    n_defining_snps += count_define_called_snps(hg_to_snps[hg], genotypes[n], issog_id_to_pos)

                hg = get_parent_hg(group_to_parent, hg)

            # record hg score
            if n_defining_snps > 0:
                all_hg_to_score[hg_names[h]] = n_called_snps / float(n_defining_snps)
                if has_parent_calls(hg_names[h], n, hg_scores, hg_names, group_to_parent, ancestral_hg_depth):
                    strict_hg_to_score[hg_names[h]] = n_called_snps / float(n_defining_snps)

    # return hg scores
    if list(strict_hg_to_score.keys()):
        return strict_hg_to_score
    else:
        return all_hg_to_score


def pick_leaf(hg_to_score, group_to_parent, outfile, sample_id, min_hap_score, min_deep_score, hg_to_snps, trunc_haps):
    """
    of the non-zero scored haplogroups collect all of the leaves, ie those which are not ancestral to any other
    non-zero haplogroup. Then, choose the group with the highest score and longest name
    """
    
    candidates = list(hg_to_score.keys())
    if not candidates:  # no hg matches, likely a poor quality sample
        print(('No match: ' + sample_id))
        outfile.write(sample_id + '\tno match\n')
        return

    leaves = set()											# leaves is now a set instead of list
    for candidate in candidates:
    	c_ancestors = get_ancestry(candidate, group_to_parent)
    	if hg_to_score[candidate] >= min_hap_score and c_ancestors:     # check if hg score is high enough and is not the root
    		is_leaf = True
    		bad_leaves = []
    		for leaf in leaves:
    			l_ancestors = get_ancestry (leaf, group_to_parent)
    			if candidate in l_ancestors:
    				is_leaf = False
    			if leaf in c_ancestors:
    				bad_leaves.append(leaf)
    		for bad_leaf in bad_leaves:
    			leaves.discard(bad_leaf)
    		if is_leaf:
    			leaves.add(candidate)
            
    max_leaf = 'A0-T'
    max_score = 0
    try:
    	max_score = hg_to_score[max_leaf]
    except:
    	max_score = 0
    if leaves: # find leaf with highest hg score
    	# just initializing variable
    	max_leaf = list(leaves)[0]
    	max_score = hg_to_score[max_leaf]
    	for leaf in leaves:
    		if max_score < hg_to_score[leaf]:
    			max_leaf = leaf
    			max_score = hg_to_score[leaf]
    	# check if there is a more derived leaf with a high score
    	longest_leaf = max(leaves, key=len)
    	if len(max_leaf) < len(longest_leaf) and hg_to_score[longest_leaf] >= min_deep_score:
    		max_leaf = longest_leaf
    		max_score = hg_to_score[longest_leaf]
    # might be better to issue a warning, then just make assignment to highest score that is most derived, or create option to just assign as root	
    elif max_score < min_hap_score:
    	print('%s: No supported leaf haplogroup available. Assigning default root haplogroup A0-T. See .all file for best assignments.' % (sample_id))    
    assign_hg = trunc_haps[max_leaf]    
    # write to output file
    hg_snps = ','.join(hg_to_snps[max_leaf])		#still showing markers for haplogroup, not the trunated haplogroup
    outfile.write('%s\t%s\t%s\t%s\n' % (sample_id, assign_hg, str(round(max_score, 3)), hg_snps))


def get_all_subgroups(hg_to_score, fi, sample_id):
    """record all non-zero scored haplogroups as a reference"""
    top_candidates = []
    candidate_scores = []
    all_scores = list(hg_to_score.values())
    while all_scores:       # order hgs based on score
        max_score = max(all_scores)

        candidates = []
        for hg in hg_to_score:
            if hg_to_score[hg] == max_score:
                candidates.append(hg)

        candidates.sort(key=len)
        candidates.reverse()

        for candidate in candidates:
            top_candidates.append(candidate)
            candidate_scores.append(max_score)

        while max_score in all_scores:
            all_scores.remove(max_score)

    # write hgs and scores to output
    line = []
    for i in range(len(top_candidates)):
        line.append(top_candidates[i] + ':' + str(round(candidate_scores[i], 3)))
    fi.write(str(sample_id) + '\t' + '\t'.join(line) + '\n')

def assign_subgroups(path, group_to_parent, samples, hg_scores, hg_to_snps, genotypes, issog_id_to_pos, sample_id, min_hap_score , min_deep_score, out_prefix, ancestral_hg_depth, trunc_haps):
    """For a sample, score all the hgs based on # derived alleles, then assign hg"""
    print('\nNow finding best-supported haplogroup for each individual')
    print('Minimum considered haplogroup score = %s' % (min_hap_score))
    print('Minimum switch to deeper node score (min_deep_score) = %s' % (min_deep_score))
    # score all hgs, then use to assign hg to individual
    with open(path + '/' + out_prefix + '.out', 'w') as leaf_outfile, open(path + '/' + out_prefix + '.all', 'w') as all_outfile:
        print('\nPrinting results to .out and .all with prefix "%s"' % (out_prefix))
        for n in range(samples):
            hg_to_score = score_hgs(hg_scores, hg_to_snps, genotypes, n, issog_id_to_pos, group_to_parent, ancestral_hg_depth)
            pick_leaf(hg_to_score, group_to_parent, leaf_outfile, sample_id[n], min_hap_score , min_deep_score, hg_to_snps, trunc_haps)
            get_all_subgroups(hg_to_score, all_outfile, sample_id[n])

# was supposed to be a recursive way of adding truncated haplogroup names
# not currently working, probably not necessary to do it recursively anyway...			        
def trunc_hap(hg, trunc_haps, group_to_parent):
	print("\tTruncating %s" % (hg))
	if hg == '':
		hg='A0-T'				#If already past the root, just set hg to root.
	if hg in ['A0-T', 'A00']:
		trunc_haps[hg] = hg
	while hg not in trunc_haps:
		parent_hg = get_parent_hg(group_to_parent, hg)
		if parent_hg == '':
			parent_hg='A0-T'
			print("Error!")
			#sys.exit()
		print("Not already stored. Now checking for %s" % (parent_hg))
		if parent_hg in trunc_haps:
			print("%s is stored. Truncating %s to %s" % (parent_hg, hg, trunc_haps[parent_hg]))
			trunc_haps[hg] = trunc_haps[parent_hg]	
		hg = parent_hg

# could be improved with a better haplogroup dictionary- fix in re-write	           
def get_trunc_haps(haps_file, all_hgs, group_to_parent):
	"""get a list of truncated haplogroup names to use in assignments"""
	trunc_haps = dict()
	if haps_file:
		with open(haps_file, 'r') as infp:
			for line in infp:
				hap = line.strip()
				trunc_haps[hap] = hap
				for parent_hap in get_ancestry(hap, group_to_parent):		#add parents of truncated haplogroups too
					if parent_hap not in trunc_haps:
						trunc_haps[parent_hap] = parent_hap
		for hg in all_hgs:		#go through each haplogroup and find its nearest ancestor from the list of truncated haplogroups
			if hg not in trunc_haps:
				ancestor_hgs = get_ancestry(hg, group_to_parent)
				for parent_hg in ancestor_hgs:
					if parent_hg in trunc_haps:
						trunc_haps[hg] = trunc_haps[parent_hg]
						break
				if not hg in trunc_haps:
					trunc_haps[hg] = hg
			#trunc_hap(hg, trunc_haps, group_to_parent)		#instead of the above loop, would try this for each hg in all_hgs if function worked...
	else:
		for hg in all_hgs:
			trunc_haps[hg] = hg
	return trunc_haps
		

def snappy(args):
	path = os.getcwd()
	project_name = args.infile
	file_prefix = path + '/' + project_name
	bed = file_prefix + '.bed'
	vcf = file_prefix + '.vcf'
	raw = file_prefix + '.raw'

	# create .raw plink file to interpret genotypes from binary files
	if not os.path.isfile(raw):
		if os.path.isfile(bed):
			print('Using plink to create .raw file from %s plink library' % (project_name))
			subprocess.call(['plink', '--bfile', project_name, '--recodeAD', '--out', project_name])
		elif os.path.isfile(vcf):
			print('Using plink to create .raw file from vcf %s' % (vcf))
    		#subprocess.call(['plink', '--vcf', project_name, '--recodeAD', '--out', project_name])
			subprocess.call(['plink', '--vcf', '%s.vcf' % (project_name), '--recodeAD', '--out', project_name])
		else:
			print('Unable to find suitable genotpye files for processing. Please ensure that there is a plink library or vcf with the prefix provided (%s)' % (file_prefix))
	else:
		print('Using %s for genotype input' % (raw))

	# build reference dictionaries
	der_allele_dict = build_derived_allele_dict('%s/%s/%s' % (path, args.ref_files_dir, args.pos2allele) )
	bim_id_dict, bim_allele_dict = build_bim_id_dict( '%s/%s.bim' % (path, project_name) )
	hg_snp_dict = build_hg_snp_dict('%s/%s/%s' % (path, args.ref_files_dir, args.hg2snp) )
	issog_id_dict = build_isogg_id_dict('%s/%s/%s' % (path, args.ref_files_dir, args.id2pos) )
	
	# read in structure of tree for non-conforming haplogroup names
	group_to_parent = getHaploGroup2Parent('%s/%s/%s' % (path, args.ref_files_dir, args.tree_strct))

	# get the genotype calls for each sample
	genotypes = []
	sample_ids = []
	n_individuals = 0
	with open(raw, 'r') as raw_data:
		bim_snp_ids = raw_data.readline().rstrip('\n').split(' ')[6:]
		for line in raw_data:
			line = line.rstrip('\n').split(' ')
			sample_id = line[1]
			sample_ids.append(sample_id)

			data = line[6:]
			genotype = get_individual_gt(bim_allele_dict, bim_id_dict, bim_snp_ids, data)
			genotypes.append(genotype)
			n_individuals += 1

	# use genotype calls to track number of derived snps called for a hg
	haplogroup_score, hg_snp_dict = tally_defining_snps(n_individuals, genotypes, hg_snp_dict, issog_id_dict, der_allele_dict)

	trunc_haps = get_trunc_haps(args.truncate_haps, list(hg_snp_dict.keys()), group_to_parent)
	# assign samples to hg
	assign_subgroups(path, group_to_parent, n_individuals, haplogroup_score, hg_snp_dict, genotypes, issog_id_dict, sample_ids, args.min_hap_score , args.min_deep_score, args.out, args.ancestral_hg_depth, trunc_haps)    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='SNAPPY')
    
    parser.add_argument('--version', action='version', version='%(prog)s 0.2.2')
    parser.add_argument('--infile', help='prefix to plink library', required=True)
    parser.add_argument('--min_hap_score', help='minimum haplogroup score to be considered as assignment', nargs='?', const=1, type=float, default=0.75, required=False)
    parser.add_argument('--min_deep_score', help='minimum score to switch to deeper node for final assignment', nargs='?', const=1, type=float, default=0.8, required=False)
    parser.add_argument('--out', help='prefix for file output', nargs='?', const=1, type=str, default='chrY_hgs', required=False)
    parser.add_argument('--ref_files_dir', help='directory where reference file are stored', nargs='?', const=1, type=str, default='ref_files', required=False)
    parser.add_argument('--id2pos', help='file listing SNP ids and corresponding positions', nargs='?', const=1, type=str, default='id_to_pos.txt', required=False)
    parser.add_argument('--pos2allele', help='file listing SNP positions and corresponding alleles', nargs='?', const=1, type=str, default='pos_to_allele.txt', required=False)
    parser.add_argument('--hg2snp', help='file listing markers and haplogroups', nargs='?', const=1, type=str, default='y_hg_and_snps.sort', required=False)
    parser.add_argument('--tree_strct', help='file listing haplogroup parent-child relationships for haplogroups that do not confrom to naming convetions', nargs='?', const=1, type=str, default='tree_structure.txt', required=False)
    parser.add_argument('--ancestral_hg_depth', help='number of ancestral haplogroups to check when considering whether a haplogroup receives a score', nargs='?', const=1, type=int, default=2, required=False)
    parser.add_argument('--truncate_haps', help='file with list of haplogroups past which SNAPPY will not make assignments', nargs='?', const=1, type=str, required=False)
    #parser.add_argument('--truncate_haps', help='file with list of haplogroups past which SNAPPY will not make assignments', action="store_const", const=1, required=False)
    
    args = parser.parse_args()
    main(args)
    sys.exit()
