import argparse
from snappy import *

def run_snappy():
        parser = argparse.ArgumentParser(prog='SNAPPY', description="Y-chromosome haplogroup inference")

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
        args = parser.parse_args()
        snappy(args)
        
if __name__ == "__main__":
        run_snappy()
