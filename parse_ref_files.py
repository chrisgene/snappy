import re


def build_bim_id_dict(filename):
    """creates two dictionaries, one to map snp ids to positions and one for positions to observed alleles"""
    id_to_pos = dict()
    pos_to_alleles = dict()
    with open(filename, 'r') as bim:
        for line in bim:
            chrom, bim_id, cm_pos, pos, allele_1, allele_2 = line.rstrip('\n').split('\t')
            if allele_1 != '0' or allele_2 != '0':
                id_to_pos[bim_id] = pos
                pos_to_alleles[pos] = allele_1 + ',' + allele_2

    return id_to_pos, pos_to_alleles


def build_hg_snp_dict(filename):
    """creates a dictionary that maps haplogroups to defining snps"""
    hg_to_snps = dict()
    with open(filename, 'r') as y_hg:
    	print('Opened %s to find haplogroup-informative SNPs.' % (filename))
    	y_hg.readline()
    	for line in y_hg:
    		hg, snps = line.rstrip('\r\n').split('\t')
    		hg_to_snps[hg] = snps.split(',')
    return hg_to_snps


def build_isogg_id_dict(filename):
    """creates a dictionary which maps isogg ids to positions from three reference files"""
    snp_id_to_pos = dict()
    with open(filename, 'r') as infp:
    	print('Creating id to position map from %s' % filename)
    	header = infp.readline()
    	for line in infp:
    		(id, pos) = line.strip('\r\n').split('\t')
    		pos = str(pos)
    		if not id in snp_id_to_pos:
    			snp_id_to_pos[id] = pos
    		else:
    			print('Warning: id %s is duplicated' % (id))

    return snp_id_to_pos


def build_derived_allele_dict(filename):
    """creates a dictionary that maps a position to a string containing the derived and ancestral alleles"""
    derived_dict = dict()
    with open(filename, 'r') as infp:
    	print('\nCreating position to allele map from %s' % filename)
    	header = infp.readline()
    	for line in infp:
    		(pos, anc, der) = line.rstrip('\r\n').split('\t')
    		pos = str(pos)
    		if not pos in derived_dict:
    			derived_dict[pos] = der + anc
    		else:
    			print('Warning: pos %s is duplicated' % (pos))
    			
    return derived_dict

def getHaploGroup2Parent(filename):
	'''create a dictionary that sends a hg to its parent hg for haplogroup names that don't conform to naming conventions'''
	group_to_parent = dict()
	with open(filename, 'r') as weird_hgs:
		for line in weird_hgs:
			line = line.rstrip('\n').split('\t')
			group_to_parent[line[0]] = line[1]
	return group_to_parent

    
