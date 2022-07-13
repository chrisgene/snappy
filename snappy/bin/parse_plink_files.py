import numpy as np


def get_individual_gt(pos_to_alleles, bim_id_to_pos, bim_snp_id, raw_gt):
    """creates a dictionary for the sample which maps position to observed genotype"""
    genotypes = dict()

    for i in range(len(raw_gt)):
        gt = raw_gt[i]
        if bim_snp_id[i][:-2] in bim_id_to_pos:     # make sure bim id in id_to_pos dic (otherwise is a nocall)
        	pos = bim_id_to_pos[bim_snp_id[i][:-2]]
        	derived_allele, ancestral_allele = pos_to_alleles[pos].split(',')
        	bim_allele = bim_snp_id[i][-1]
        	
        	# record genotype call in dictionary
        	if bim_allele == derived_allele:
        		if gt == '0':
        			genotypes[pos] = ancestral_allele
        		elif gt == '2':
        			genotypes[pos] = derived_allele
        		else:
        			genotypes[pos] = 'N'
        	elif bim_allele == ancestral_allele:
        		if gt == '0':
        			genotypes[pos] = derived_allele
        		elif gt == '2':
        			genotypes[pos] = ancestral_allele
        		else:
        			genotypes[pos] = 'N'
        	elif bim_allele == '0':
        		if gt == '0':
        			genotypes[pos] = ancestral_allele
        
    return genotypes


def tally_defining_snps(n_samples, genotypes, hg_to_snp, id_to_pos, pos_to_derived_allele):
    """for each haplogroup, count how many of the defining snps are derived and record efficiently in numpy array"""
    hg_score = np.zeros((n_samples, len(hg_to_snp)))
    hg_names = list(hg_to_snp.keys())

    for n in range(n_samples):
        pos_to_gt = genotypes[n]
        for j in range(len(hg_names)):
            defining_snps = hg_to_snp[hg_names[j]]
            n_defining_snps = 0
            for snp in defining_snps:
                # check if snp has multiple ids, find one which will give position
                if '/' in snp:
                    same_snp = snp.split('/')
                    for s in same_snp:
                        if s in id_to_pos:
                            snp = s
                            break
                # count snp if the allele is derived
                if snp in id_to_pos:
                    pos = id_to_pos[snp]
                    if pos in pos_to_gt and pos in pos_to_derived_allele:
                    #if pos in pos_to_gt:
                        gt = pos_to_gt[pos]
                        derived_allele = pos_to_derived_allele[pos][0]
                        if gt == derived_allele:
                            n_defining_snps += 1
                    else:
                        # position not in bim file, need if statement bc will try to delete multiple times
                        if snp in hg_to_snp[hg_names[j]]:
                            hg_to_snp[hg_names[j]].remove(snp)

                else:
                    # don't know position of defining snp
                    hg_to_snp[hg_names[j]].remove(snp)

            hg_score[n, j] = n_defining_snps

    return hg_score, hg_to_snp
