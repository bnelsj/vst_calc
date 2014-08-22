#! /usr/bin python

import numpy as np
import pandas as pd
import argparse
import pysam
import cluster
import genotyper as gt
from GC_data import GC_data

def get_pops_from_file(fn_pops):
    pops = {}
    
    with open(args.populations, 'r') as pops_reader:
        for line in pops_reader:
            indiv, sex, pop = line.split()
            if pop not in pops:
                pops[pop] = []
            pops[pop].append(indiv)
    return(pops)

def calc_vst(pop1, pop2):
    #Calculate vst between two populations
    pop1 = pop1[pop1 != -1]
    pop2 = pop2[pop2 != -1]
    n1 = len(pop1)
    n2 = len(pop2)
    
    if n1 == 0 or n2 == 0 or (n1 + n2 < 30):
         return 0

    n_total = n1 + n2
    v_total = np.var(np.r_[pop1, pop2])
    if v_total == 0: return 0

    v1 = np.var(pop1)
    v2 = np.var(pop2)
    vst = (v_total - (((v1*n1)+(v2*n2))/n_total))/v_total
    return vst

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gglob_dir", default="/net/eichler/vol22/projects/1000_genomes_phase_II_III/nobackups/jlhudd/dCGH/2014-02-05-1kg_vs_HGDP_redux/gglob", 
                        help="Path to gglob directory (Default: %(default)s)")
    parser.add_argument("--populations", default="/net/eichler/vol22/projects/1000_genomes_phase_II_III/nobackups/bnelsj/1kg_genotypes/vst_stats/sex.genome.1KG_PHASE_I_II_III_GENERAL", 
                        help="Path to file with sample->population info (Default: %(default)s)")   
    parser.add_argument("--plot_dir", default="plots", help="Path to plotting directory")
    parser.add_argument("--fn_fa", default="/net/eichler/vol7/home/psudmant/genomes/fastas/hg19_1kg_phase2_reference/human_g1k_v37.fasta")
    #"/net/eichler/vol7/home/psudmant/genomes/annotations/hg19/superdups/superdups.merged.bed.gz"
    parser.add_argument("--fn_dup_tabix", default="/net/eichler/vol7/home/psudmant/genomes/annotations/hg19/superdups/superdups.merged.bed.gz")
    parser.add_argument("--fn_GC_DTS", default="/net/eichler/vol7/home/psudmant/genomes/GC_tracks/windowed_DTS/HG19/500_bp_slide_GC")
    parser.add_argument("--fn_DTS_contigs", default="/net/eichler/vol7/home/psudmant/EEE_Lab/1000G/1000genomesScripts/windowed_analysis/DTS_window_analysis/windows/hg19_slide/500_bp_windows.pkl.contigs")
    parser.add_argument("--regions", default="/net/eichler/vol2/eee_shared/assemblies/hg19/genes/refGene.bed", help="Path to regions bed file (Default: %(default)s)")
    parser.add_argument("--contig", required=True, help="Chromosome to consider")
    parser.add_argument("--outfile", required=True, help="Path to output file")
    parser.add_argument("--header", help="Chromosome header will be printed for")

    args = parser.parse_args()

    pop_dict = sorted(get_pops_from_file(args.populations))

    pops = pd.read_csv(args.populations, header=None, delimiter="\t", index_col=None)
    pops.columns = ["indiv", "sex", "pop"]
    pops_unique = pd.unique(pops["pop"])
    pop_combos = [(x, y) for i, x in enumerate(pops_unique) for j, y in enumerate(pops_unique) if i < j]
    pop_header = map(lambda x: "-".join(x), pop_combos)

    indivs = list(pd.read_json("%s/gglob.idx" % args.gglob_dir).indivs)

    pop_indices = {}

    for i, indiv in enumerate(indivs):
        pop = pops[pops["indiv"] == indiv]["pop"].values[0]
        if pop not in pop_indices:
            pop_indices[pop] = []
        pop_indices[pop].append(i)

    regions = pd.read_csv(args.regions, header=None, delimiter="\t", index_col=None)
    ncols_regions = regions.shape[1]
    if ncols_regions == 4:
        regions.columns = ["contig", "start", "end", "name"]
    elif ncols_regions == 6:
        regions.columns = ["contig", "start", "end", "name", "score", "strand"]
    else:
        print "Regions file should have 4 or 6 columns, not %d" % ncols_regions

    regions_subset = regions[regions["contig"] == args.contig]

    vst_cols = ["contig", "start", "end", "name"].extend(pop_combos)

    tbx_dups = pysam.Tabixfile(args.fn_dup_tabix)
    GC_inf = GC_data(args.fn_GC_DTS, args.contig, args.fn_DTS_contigs)

    g = gt.genotyper(args.contig, gglob_dir = args.gglob_dir, plot_dir = args.plot_dir, subset_indivs = indivs, fn_fa=args.fn_fa, dup_tabix = tbx_dups, GC_inf = GC_inf)

    FOUT = open(args.outfile, 'w')

    if args.header is not None and args.header == args.contig:
        print "printing header...\n"
        FOUT.write("contig\tstart\tend\tname\t%s\n"%("\t".join(pop_header)))

    for i, row in regions_subset.iterrows():
        contig, s, e, name = row['contig'], row['start'], row['end'], row['name']
        X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        cps = np.mean(X,1)
        outrow = [contig, s, e, name]
        for pop_pair in pop_combos:
            cps0 = cps[pop_indices[pop_pair[0]]]
            cps1 = cps[pop_indices[pop_pair[1]]]
            vst = calc_vst(cps0, cps1)
            outrow.append(vst)
        FOUT.write("\t".join(map(str,outrow)) + "\n")

    FOUT.close()
