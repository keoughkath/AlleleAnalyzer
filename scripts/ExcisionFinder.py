# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ExcisionFinder identifies allele-specific excision sites. Written in Python version 3.6.1.
Kathleen Keough et al 2017-2018.

Usage: 
        ExcisionFinder.py [-vsgc] <gene_dat> <gene> <annots_file> <maxcut> <cas_list> <bcf> <out> [--window=<window_in_bp>] [--not_phased] [--guides=<guides>] [--exhaustive]
        ExcisionFinder.py -h

Arguments:
    gene_dat                         Gene annotations file (gene_gene_dat_wsize) filepath.
    gene                             Gene you would like to analyze.
    var_annots                       Variant annotation HDF5 file.
    maxcut                           Maximum distance between cut position pairs.
    cas_list                         Comma separated (no spaces!) list of Cas varieties to evaluate, options below.
    bcf                              BCF/VCF file with variants for individual or cohort being analyzed. 
    out                              Directory to which you would like to write the output files.

Options:
    -h                               Show this screen.
    -v                               Run as verbose, print to stderr. 
    -s                               Only consider sgRNA sites where variant is in a PAM (strict).
    --window=<window_in_bp>          Window around the gene (in basepairs) to also consider [default: 0]. 
    -g                               Output guides for heterozygous variants that generate targetable pairs. Requires <guides> file.
    --not_phased                     Add this if genomes are not phased to skip evaluating whether sites are
                                     on same haplotype.
    --guides=<guides>                Guides file for locus if '-g' specified.
    --exhaustive                     Run exhaustive style analysis (e.g. for set cover analysis)

Available Cas types = cpf1,SpCas9,SpCas9_VRER,SpCas9_EQR,SpCas9_VQR_1,SpCas9_VQR_2,StCas9,StCas9_2,SaCas9,SaCas9_KKH,nmCas9,cjCas9
"""

import pandas as pd
from pandas import HDFStore
import numpy as np
from functools import reduce
from docopt import docopt
import itertools
import regex as re
import logging
import subprocess
from io import StringIO
import os, sys
import time
import cas_object as cas_obj

# Get absolute path for ExcisionFinder.py, and edit it for cas_object.py
ef_path = os.path.dirname(os.path.realpath(__file__))
metadata_path = '/'.join(ef_path.split('/')[:-1] + ['preprocessing'])
sys.path.append(metadata_path)
from get_metadata import add_metadata

__version__ = '0.0.0'


def load_gene_gene_dat(gene_dat_path):
    """
    Load gene annotation data (transcript data).
    :param gene_dat_path: str, filepath for gene_gene_dat_wsize (Part of ExcisionFinder package).
    :return: Refseq gene annotations file.
    """
    gene_gene_dat = pd.read_csv(
        gene_dat_path,
        sep="\t",
        header=0,
        names=[
            "name",
            "chrom",
            "txStart",
            "txEnd",
            "cdsStart",
            "cdsEnd",
            "exonCount",
            "exonStarts",
            "exonEnds",
            "size",
        ],
    )
    return gene_gene_dat


def het(genotype):
    """
    Determine whether a genotype in format A|G is het.
    :param genotype: genotype, str.
    :return: bool, True = het, False = hom.
    """
    if ':' in genotype:
        genotype = genotype.split(':')[0]
    hap1, hap2 = re.split("/|\|", genotype)
    return hap1 != hap2


def next_exon(variant_position, coding_exon_starts):
    """
    get location of next coding exon after variant
    :param coding_exon_starts: coding exon start positions, Pandas Series.
    :param variant_position: chromosomal position, int
    :return: chromosomal position of start of next coding exon, int
    """
    greater_than_var = [x for x in coding_exon_starts if x > variant_position]
    if not greater_than_var:
        return False
    else:
        next_coding_exon_pos = min(greater_than_var)
        return next_coding_exon_pos


def targ_pair(variant1, variant2, coding_positions, coding_exon_starts):
    """
    Determine whether a pair of variants positions is targetable based on whether they might
    disrupt an exon.
    :param variant1: position of variant 1, int.
    :param variant2: position of variant 2, int.
    :param coding_positions: coding positions, set.
    :param coding_exon_starts: Start positions of coding exons, Pandas Series.
    :return: whether targetable or not, bool.
    """
    low_var, high_var = sorted([variant1, variant2])
    if low_var in coding_positions or high_var in coding_positions:
        return True
    elif next_exon(low_var, coding_exon_starts):
        # checks whether larger variant position occurs in or after next exon
        return bool(high_var >= next_exon(low_var, coding_exon_starts))
    else:
        return False


def translate_gene_name(gene_name):
    """
    HDF5 throws all sort of errors when you have weird punctuation in the gene name, so
    this translates it to a less offensive form.
    """
    repls = ("-", "dash"), (".", "period")
    trans_gene_name = reduce(lambda a, kv: a.replace(*kv), repls, str(gene_name))
    return trans_gene_name


class Gene:
    """Holds information for the gene"""

    def __init__(self, official_gene_symbol, gene_dat, window):
        self.official_gene_symbol = official_gene_symbol
        self.info = gene_dat.query("index == @self.official_gene_symbol")
        self.n_exons = self.info["exonCount"].item()
        self.coding_start = int(self.info["cdsStart"].item())
        self.coding_end = int(self.info["cdsEnd"].item())
        self.coding_exons = []
        counter = 0
        counterweight = len(list(
                zip(
                    list(map(int, self.info["exonStarts"].item().split(",")[:-1])),
                    list(map(int, self.info["exonEnds"].item().split(",")[:-1])),
                )))
        for x in list(zip(list(map(int, self.info["exonStarts"].item().split(",")[:-1])),
                      list(map(int, self.info["exonEnds"].item().split(",")[:-1])))):
            counter += 1
            if counter == 1:
                self.coding_exons.append((max(self.coding_start, x[0]), x[1]))
            elif counter == counterweight:
                self.coding_exons.append((x[0], min(self.coding_end, x[1])))
            else:
                self.coding_exons.append((x[0], x[1]))
        self.n_coding_exons = len(self.coding_exons)
        self.start = self.info["txStart"].item() - window
        self.end = self.info["txEnd"].item() + window
        self.chrom = gene_dat.query("index == @self.official_gene_symbol")[
            "chrom"
        ].item()

    def get_coding_positions_and_starts(self):
        coding_positions = []
        coding_exon_starts = []
        for start, stop in self.coding_exons:
            coding_positions.extend(list(range(start, stop + 1)))
            coding_exon_starts.append(start)
        return coding_positions, coding_exon_starts


def check_bcftools():
    """ 
    Checks bcftools version, and exits the program if the version is incorrect
    """
    version = (
        subprocess.run(
            "bcftools -v | head -1 | cut -d ' ' -f2", shell=True, stdout=subprocess.PIPE
        )
        .stdout.decode("utf-8")
        .rstrip()
    )
    if float(version) >= REQUIRED_BCFTOOLS_VER:
        print(f"bcftools version {version} running")

    else:
        print(f"Error: bcftools must be >=1.5. Current version: {version}")
        exit()


def pair_guides(guides_df, variant1, variant2):
    guides_df = pd.read_csv(guides_df, sep="\t")
    guide_pairs_out = pd.DataFrame()
    guide_pairs_out["variant1"] = variant1
    guide_pairs_out["variant2"] = variant2
    guide_pairs_out["variant_position"] = guide_pairs_out["variant1"]
    guide_pairs_out = guide_pairs_out.merge(
        guides_df, how="left", on=["variant_position"]
    )
    guide_pairs_out["gRNA_ref_guide_1"] = guide_pairs_out["gRNA_ref"]
    guide_pairs_out["gRNA_alt_guide_1"] = guide_pairs_out["gRNA_alt"]
    guide_pairs_out["ref_g1"] = guide_pairs_out["ref"]
    guide_pairs_out["alt_g1"] = guide_pairs_out["alt"]
    guide_pairs_out['id_g1'] = guide_pairs_out['id']
    guide_pairs_out["variant_position_in_guide_g1"] = guide_pairs_out[
        "variant_position_in_guide"
    ]
    guide_pairs_out = guide_pairs_out[
        [
            "variant1",
            "variant2",
            "gRNA_ref_guide_1",
            "gRNA_alt_guide_1",
            "ref_g1",
            "alt_g1",
            "variant_position_in_guide_g1",
            "id_g1"
        ]
    ]
    guide_pairs_out["variant_position"] = guide_pairs_out["variant2"]
    guide_pairs_out = guide_pairs_out.merge(
        guides_df, how="left", on=["variant_position"]
    )
    guide_pairs_out["gRNA_ref_guide_2"] = guide_pairs_out["gRNA_ref"]
    guide_pairs_out["gRNA_alt_guide_2"] = guide_pairs_out["gRNA_alt"]
    guide_pairs_out["ref_g2"] = guide_pairs_out["ref"]
    guide_pairs_out["alt_g2"] = guide_pairs_out["alt"]
    guide_pairs_out['id_g2'] = guide_pairs_out['id']
    guide_pairs_out["variant_position_in_guide_g2"] = guide_pairs_out[
        "variant_position_in_guide"
    ]
    guides_out = guide_pairs_out[
        [
            "variant1",
            "variant2",
            "gRNA_ref_guide_1",
            "gRNA_alt_guide_1",
            "ref_g1",
            "alt_g1",
            "variant_position_in_guide_g1",
            "id_g1",
            "gRNA_ref_guide_2",
            "gRNA_alt_guide_2",
            "ref_g2",
            "alt_g2",
            "variant_position_in_guide_g2",
            "id_g2"
        ]
    ].dropna()
    return guides_out


def targ_var(row, cas, hap, level="lax"):
    var = row["pos"]
    if level == "strict":
        if hap == 1:
            return annots_file.query("pos == @var")[[f"makes_{cas}"]].values.any()
        elif hap == 2:
            return annots_file.query("pos == @var")[[f"breaks_{cas}"]].values.any()
    else:
        if hap == 1:
            return annots_file.query("pos == @var")[
                [f"makes_{cas}", f"var_near_{cas}"]
            ].values.any()
        if hap == 2:
            return annots_file.query("pos == @var")[
                [f"breaks_{cas}", f"var_near_{cas}"]
            ].values.any()


def norm_chr(chrom_str, gens_chrom):
    """
    Returns the chromosome string that matches the chromosome annotation of the input gens file
    """
    chrom_str = str(chrom_str)
    if not gens_chrom:
        return chrom_str.replace("chr", "")
    elif gens_chrom and not chrom_str.startswith("chr"):
        return "chr" + chrom_str
    else:
        return chrom_str


def main(args):

    gene_dat = load_gene_gene_dat(args["<gene_dat>"])
    gene = args["<gene>"]
    global annots_file
    annots_file = args["<annots_file>"]
    out_prefix = args["<out>"]
    maxcut = int(args["<maxcut>"])
    cas_list_append = args["<cas_list>"].split(",")
    bcf = args["<bcf>"]
    window = int(args["--window"])

    cas_list = ["all"] + cas_list_append

    # define strictness level, which is whether or not variants near PAMs are considered
    # along with those that are in PAMs

    if args["-s"]:
        logging.info("Running as strict.")
        strict_level = "strict"
    else:
        strict_level = "relaxed"
        logging.info("Running as relaxed.")

    logging.info("Now running ExcisionFinder on " + gene + ".")

    # grab info about relevant gene w/ class

    MyGene = Gene(gene, gene_dat, window)

    # get number of coding exons in gene, must have at least 1 to continue

    n_exons = MyGene.n_exons
    n_coding_exons = MyGene.n_coding_exons

    if n_coding_exons < 1:
        logging.error(
            f"{n_exons} total exons in this gene, {n_coding_exons} of which are coding.\
            No coding exons in gene {gene}, exiting."
        )
        with open(f"{out_prefix}no_coding_exons.txt", "a+") as f:
            f.write(gene + "\n")
        exit()
    else:
        logging.info(
            f"{n_exons} total exons in this gene, {n_coding_exons} of which are coding."
        )

    # load targetability information for each variant
    annots_file = pd.read_hdf(
        annots_file, "all", where=f"pos >= {MyGene.start} and pos <= {MyGene.end}"
    )

    # check whether there are annotated variants for this gene, abort otherwise

    if annots_file.empty:
        logging.error(f"No variants in file for gene {gene}")
        with open(f"{out_prefix}not_enough_hets.txt", "a+") as fout:
            fout.write(gene + "\n")
        exit()
    else:
        logging.info(
            f"Targetability data loaded, {annots_file.shape[0]} variants annotated in 1KGP for {gene}."
        )

    # import region of interest genotypes

    bcf = f"{bcf}"

    vcf_chrom = str(
        subprocess.Popen(
            f"bcftools view -H {bcf} | cut -f1 | head -1",
            shell=True,
            stdout=subprocess.PIPE,
        ).communicate()[0].decode("utf-8")
    )
    # See if chrom contains chr

    chrom = norm_chr(MyGene.chrom, vcf_chrom.startswith("chr"))

    bcl_v = f'bcftools view -g "het" -r {chrom}:{MyGene.start}-{MyGene.end} -H {bcf}'

    samples_cmd = f"bcftools query -l {bcf}"
    bcl_samps = subprocess.Popen(samples_cmd, shell=True, stdout=subprocess.PIPE)
    samples = bcl_samps.communicate()[0].decode("utf-8").split("\n")[:-1]

    # check that user specified cohort if the VCF contains >1 samples

    col_names = [
        "chrom",
        "pos",
        "rsid",
        "ref",
        "alt",
        "score",
        "random",
        "info",
        "gt",
    ] + samples
    bcl_view = subprocess.Popen(bcl_v, shell=True, stdout=subprocess.PIPE)
    gens = pd.read_csv(
        StringIO(bcl_view.communicate()[0].decode("utf-8")),
        sep="\t",
        header=None,
        names=col_names,
        usecols=["chrom", "pos", "ref", "alt"] + samples,
    )
    logging.info("Genotypes loaded.")

    het_gens = gens[samples].applymap(het).copy()

    enough_hets = list(het_gens.sum(axis=0).loc[lambda s: s >= 2].index)

    logging.info(str(len(enough_hets)) + " individuals have >= 2 het positions.")

    if len(enough_hets) < 1:
        logging.info("No individuals have at least 2 het sites, aborting analysis.")
        with open(f"{out_prefix}not_enough_hets.txt", "a+") as fout:
            fout.write(gene + "\n")
        exit()

    logging.info(
        "Checking targetability of individuals with sufficient number of hets."
    )

    # get variants in region

    variants = sorted(gens.pos)

    # set up targetability analyses

    het_vars_per_ind = {}  # get heterozygous variant positions for each individual

    for ind in enough_hets:
        het_vars_per_ind[ind] = gens.pos[het_gens[ind]].tolist()

    # get variant combinations and extract targetable pairs

    logging.info("Getting variant combos.")

    variant1 = []
    variant2 = []

    coding_positions, coding_exon_starts = MyGene.get_coding_positions_and_starts()

    for var1, var2 in itertools.product(variants, repeat=2):
        if (
            (var1 != var2)
            and (max([var1, var2]) <= min([var1, var2]) + maxcut)
            and (targ_pair(var1, var2, coding_positions, coding_exon_starts))
        ):
            variant1.append(var1)
            variant2.append(var2)
        else:
            continue

    logging.info("Combos obtained.")

    # get rid of dups and make df

    targ_pairs_df = pd.DataFrame({"var1": variant1, "var2": variant2}).query(
        "var1 < var2"
    )

    # check that each individual that has enough hets also has at least one of these pairs
    # and specify which pairs they have

    inds_w_targ_pair = {}

    for ind in het_vars_per_ind.keys():
        ind_vars = het_vars_per_ind[ind]
        ind_targ_pairs = (
            targ_pairs_df.loc[targ_pairs_df.isin(ind_vars).all(axis=1)]
            .reset_index(drop=True)
            .copy()
        )
        if ind_targ_pairs.empty:
            continue
        else:
            inds_w_targ_pair[ind] = ind_targ_pairs

    logging.info(
        f"{len(inds_w_targ_pair.keys())} individuals have at least one targetable pair of variants."
    )

    if not inds_w_targ_pair:
        logging.info(
            f"No individuals in 1KGP have at least 1 targetable variant pair for {gene}."
        )
        with open(f"{out_prefix}no_targetable_inds.txt", "a+") as fout:
            fout.write(gene + "\n")
        exit()

    # check targetability for each type of Cas

    final_targ = pd.DataFrame({"sample": list(inds_w_targ_pair.keys())})

    finaltargcols = []  # keeps track of columns for all cas types for later evaluating "all" condition

    if args["--exhaustive"]:
        # val_occurrences = []  # steps towards longform dataframe with one row per ind/het combo
        # for val in het_vars_per_ind.values():
        #     val_occurrences.append(len(val))
        # targ_haps = pd.DataFrame({'sample': np.repeat(list(het_vars_per_ind.keys()), val_occurrences),
        #                           'pos': list(itertools.chain.from_iterable(het_vars_per_ind.values()))})
        overall_exh_list = []

    for cas in cas_list[1:]:  # skip all because is handled below faster
        logging.info(f"Evaluating gene targetability for {cas}")
        if args["-s"]:
            targ_vars_cas = annots_file.query(
                f"(makes_{cas}) or (breaks_{cas})"
            ).pos.tolist()
        else:
            targ_vars_cas = annots_file.query(
                f"(var_near_{cas}) or (makes_{cas}) or (breaks_{cas})"
            ).pos.tolist()
        targ_pairs_cas = (
            targ_pairs_df.loc[targ_pairs_df.isin(targ_vars_cas).all(axis=1)]
            .reset_index(drop=True)
            .copy()
        )
        if args["--exhaustive"]:
            exh_df_list = []
        # eliminate individuals that do not have at least one targetable pair for this specific cas
        ind_targ_cas = []
        for ind in list(inds_w_targ_pair.keys()):
            ind_cas_targ_pairs = (
                inds_w_targ_pair[ind]
                .merge(targ_pairs_cas, how="inner")
                .drop_duplicates()
            )
            if args["--exhaustive"]:
                ind_targ_out = []
                # check whether pairs of allele-specific cut sites is on the same haplotype in individual
                if args["-s"]:
                    ind_cas_targ_pairs["var1_make_pam"] = ind_cas_targ_pairs[
                        "var1"
                    ].isin(annots_file.query(f"makes_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var2_make_pam"] = ind_cas_targ_pairs[
                        ["var2"]
                    ].isin(annots_file.query(f"makes_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var1_break_pam"] = ind_cas_targ_pairs[
                        ["var1"]
                    ].isin(annots_file.query(f"breaks_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var2_break_pam"] = ind_cas_targ_pairs[
                        ["var2"]
                    ].isin(annots_file.query(f"breaks_{cas}").pos.tolist())
                    ind_cas_targ_pairs["both_make"] = ind_cas_targ_pairs[
                        ["var1_make_pam", "var2_make_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["both_break"] = ind_cas_targ_pairs[
                        ["var1_break_pam", "var2_break_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["one_make_one_break_1"] = ind_cas_targ_pairs[
                        ["var1_make_pam", "var2_break_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["one_make_one_break_2"] = ind_cas_targ_pairs[
                        ["var2_make_pam", "var1_break_pam"]
                    ].all(axis=1)
                    gens_replace = {
                        "0|1": "hap2",
                        "0|2": "hap2",
                        "0|3": "hap2",
                        "1|0": "hap1",
                        "2|0": "hap1",
                        "3|0": "hap1",
                        "0|0": "not_het",
                        "1|1": "not_het",
                    }
                    ind_gens = gens[["pos", ind]].replace(gens_replace).copy()
                    ind_cas_targ_pairs["var1_hap"] = pd.merge(
                        ind_cas_targ_pairs.copy(),
                        ind_gens.copy(),
                        left_on="var1",
                        right_on="pos",
                        how="left",
                    )[ind]
                    ind_cas_targ_pairs["var2_hap"] = pd.merge(
                        ind_cas_targ_pairs.copy(),
                        ind_gens.copy(),
                        left_on="var2",
                        right_on="pos",
                        how="left",
                    )[ind]
                    ind_cas_targ_pairs["same_hap"] = np.where(
                        ind_cas_targ_pairs["var1_hap"]
                        == ind_cas_targ_pairs["var2_hap"],
                        True,
                        False,
                    )
                    ind_cas_targ_pairs["not_same_hap"] = ~ind_cas_targ_pairs["same_hap"]

                    ind_targ_out.append(
                        ind_cas_targ_pairs.query("both_make and same_hap")[
                            ["var1", "var2"]
                        ]
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query("both_break and same_hap")[
                            ["var1", "var2"]
                        ]
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query(
                            "one_make_one_break_1 and not_same_hap"
                        )[["var1", "var2"]]
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query(
                            "one_make_one_break_2 and not_same_hap"
                        )[["var1", "var2"]]
                    )
                    ind_targ_out_df = pd.concat(ind_targ_out).dropna().drop_duplicates()
                    ind_targ_out_df["ind"] = ind
                    exh_df_list.append(ind_targ_out_df)
                else:
                    # if both near PAM, haplotype doesn't have to be the same because both are allele-specific sgRNA sites
                    ind_cas_targ_pairs["both_near_pam"] = (
                        ind_cas_targ_pairs[["var1", "var2"]]
                        .isin(annots_file.query(f"var_near_{cas}").pos.tolist())
                        .all(axis=1)
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query("both_near_pam")[["var1", "var2"]]
                    )
                    #   when both make or break, need to be same hap
                    ind_cas_targ_pairs["var1_make_pam"] = ind_cas_targ_pairs[
                        "var1"
                    ].isin(annots_file.query(f"makes_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var2_make_pam"] = ind_cas_targ_pairs[
                        "var2"
                    ].isin(annots_file.query(f"makes_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var1_near_pam"] = ind_cas_targ_pairs[
                        "var1"
                    ].isin(annots_file.query(f"var_near_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var2_near_pam"] = ind_cas_targ_pairs[
                        "var2"
                    ].isin(annots_file.query(f"var_near_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var1_break_pam"] = ind_cas_targ_pairs[
                        "var1"
                    ].isin(annots_file.query(f"breaks_{cas}").pos.tolist())
                    ind_cas_targ_pairs["var2_break_pam"] = ind_cas_targ_pairs[
                        "var2"
                    ].isin(annots_file.query(f"breaks_{cas}").pos.tolist())
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query(
                            "(var1_near_pam and var2_make_pam) or (var1_near_pam and var2_break_pam) or (var2_near_pam and var1_make_pam) or (var2_near_pam and var1_break_pam)"
                        )[["var1", "var2"]]
                    )
                    ind_cas_targ_pairs["both_make"] = ind_cas_targ_pairs[
                        ["var1_make_pam", "var2_make_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["both_break"] = ind_cas_targ_pairs[
                        ["var1_break_pam", "var2_break_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["one_make_one_break_1"] = ind_cas_targ_pairs[
                        ["var1_make_pam", "var2_break_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["one_make_one_break_2"] = ind_cas_targ_pairs[
                        ["var2_make_pam", "var1_break_pam"]
                    ].all(axis=1)
                    gens_replace = {
                        "0|1": "hap2",
                        "0|2": "hap2",
                        "0|3": "hap2",
                        "1|0": "hap1",
                        "2|0": "hap1",
                        "3|0": "hap1",
                        "0|0": "not_het",
                        "1|1": "not_het",
                    }
                    ind_gens = gens[["pos", ind]].replace(gens_replace).copy()
                    ind_cas_targ_pairs["var1_hap"] = pd.merge(
                        ind_cas_targ_pairs.copy(),
                        ind_gens.copy(),
                        left_on="var1",
                        right_on="pos",
                        how="left",
                    )[ind]
                    ind_cas_targ_pairs["var2_hap"] = pd.merge(
                        ind_cas_targ_pairs.copy(),
                        ind_gens.copy(),
                        left_on="var2",
                        right_on="pos",
                        how="left",
                    )[ind]
                    ind_cas_targ_pairs["same_hap"] = np.where(
                        ind_cas_targ_pairs["var1_hap"]
                        == ind_cas_targ_pairs["var2_hap"],
                        True,
                        False,
                    )
                    ind_cas_targ_pairs["not_same_hap"] = ~ind_cas_targ_pairs["same_hap"]
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query("both_make and same_hap")[
                            ["var1", "var2"]
                        ]
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query("both_break and same_hap")[
                            ["var1", "var2"]
                        ]
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query(
                            "one_make_one_break_1 and not_same_hap"
                        )[["var1", "var2"]]
                    )
                    ind_targ_out.append(
                        ind_cas_targ_pairs.query(
                            "one_make_one_break_2 and not_same_hap"
                        )[["var1", "var2"]]
                    )
                    ind_targ_out_df = pd.concat(ind_targ_out).dropna().drop_duplicates()
                    ind_targ_out_df["ind"] = ind
                    exh_df_list.append(ind_targ_out_df)
                exh_df = pd.concat(exh_df_list)
                exh_df[f"targ_{cas}"] = cas
                overall_exh_list.append(exh_df)
            if ind_cas_targ_pairs.empty:
                ind_targ_cas.append(False)
                continue
            # don't need to check that pairs are on same haplotype if genotypes are not phased, unless if outputting hap_targs
            elif args["--not_phased"] and not args["--exhaustive"]:
                ind_targ_cas.append(True)
                continue
            else:
                # check that at least one pair of allele-specific cut sites is on the same haplotype in individual
                if args["-s"]:
                    ind_cas_targ_pairs["both_make"] = ind_cas_targ_pairs[
                        ["var1_make_pam", "var2_make_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["both_break"] = ind_cas_targ_pairs[
                        ["var1_break_pam", "var2_break_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["one_make_one_break_1"] = ind_cas_targ_pairs[
                        ["var1_make_pam", "var2_break_pam"]
                    ].all(axis=1)
                    ind_cas_targ_pairs["one_make_one_break_2"] = ind_cas_targ_pairs[
                        ["var2_make_pam", "var1_break_pam"]
                    ].all(axis=1)
                    gens_replace = {
                        "0|1": "hap2",
                        "0|2": "hap2",
                        "0|3": "hap2",
                        "1|0": "hap1",
                        "2|0": "hap1",
                        "3|0": "hap1",
                        "0|0": "not_het",
                        "1|1": "not_het",
                    }
                    ind_gens = gens[["pos", ind]].replace(gens_replace).copy()
                    ind_cas_targ_pairs["var1_hap"] = pd.merge(
                        ind_cas_targ_pairs.copy(),
                        ind_gens.copy(),
                        left_on="var1",
                        right_on="pos",
                        how="left",
                    )[ind]
                    ind_cas_targ_pairs["var2_hap"] = pd.merge(
                        ind_cas_targ_pairs.copy(),
                        ind_gens.copy(),
                        left_on="var2",
                        right_on="pos",
                        how="left",
                    )[ind]
                    ind_cas_targ_pairs["same_hap"] = np.where(
                        ind_cas_targ_pairs["var1_hap"]
                        == ind_cas_targ_pairs["var2_hap"],
                        True,
                        False,
                    )
                    ind_cas_targ_pairs["not_same_hap"] = ~ind_cas_targ_pairs["same_hap"]
                    if (
                        ind_cas_targ_pairs[["both_make", "same_hap"]].all(axis=1).any()
                        or ind_cas_targ_pairs[["both_break", "same_hap"]]
                        .all(axis=1)
                        .any()
                    ):
                        ind_targ_cas.append(True)
                        continue
                    # check if pair where one makes, one breaks a PAM, and on different haplotypes
                    elif (
                        ind_cas_targ_pairs[["one_make_one_break_1", "not_same_hap"]]
                        .all(axis=1)
                        .any()
                        or ind_cas_targ_pairs[["one_make_one_break_2", "not_same_hap"]]
                        .all(axis=1)
                        .any()
                    ):
                        ind_targ_cas.append(True)
                        continue
                    # all possibilities exhausted, this person just isn't targetable at this gene
                    else:
                        ind_targ_cas.append(False)
                        continue
                else:
                    # if both near PAM, haplotype doesn't have to be the same because both are allele-specific sgRNA sites
                    ind_cas_targ_pairs["both_near_pam"] = (
                        ind_cas_targ_pairs[["var1", "var2"]]
                        .isin(annots_file.query(f"var_near_{cas}").pos.tolist())
                        .all(axis=1)
                    )
                    if ind_cas_targ_pairs[
                        "both_near_pam"
                    ].any():  # this doesn't work for "strict" mode
                        ind_targ_cas.append(True)
                        continue
                    else:
                        # if none have both near a PAM, when both make or break, need to be same hap
                        ind_cas_targ_pairs["var1_make_pam"] = ind_cas_targ_pairs[
                            "var1"
                        ].isin(annots_file.query(f"makes_{cas}").pos.tolist())
                        ind_cas_targ_pairs["var2_make_pam"] = ind_cas_targ_pairs[
                            ["var2"]
                        ].isin(annots_file.query(f"makes_{cas}").pos.tolist())
                        ind_cas_targ_pairs["var1_near_pam"] = ind_cas_targ_pairs[
                            "var1"
                        ].isin(annots_file.query(f"var_near_{cas}").pos.tolist())
                        ind_cas_targ_pairs["var2_near_pam"] = ind_cas_targ_pairs[
                            "var2"
                        ].isin(annots_file.query(f"var_near_{cas}").pos.tolist())
                        ind_cas_targ_pairs["var1_break_pam"] = ind_cas_targ_pairs[
                            ["var1"]
                        ].isin(annots_file.query(f"breaks_{cas}").pos.tolist())
                        ind_cas_targ_pairs["var2_break_pam"] = ind_cas_targ_pairs[
                            ["var2"]
                        ].isin(annots_file.query(f"breaks_{cas}").pos.tolist())
                        # if one var is near a pam and the other makes/breaks, haplotype doesn't matter
                        if not ind_cas_targ_pairs.query(
                            "(var1_near_pam and var2_make_pam) or (var1_near_pam and var2_break_pam) or (var2_near_pam and var1_make_pam) or (var2_near_pam and var1_break_pam)"
                        ).empty:
                            ind_targ_cas.append(True)
                            continue
                        else:
                            ind_cas_targ_pairs["both_make"] = ind_cas_targ_pairs[
                                ["var1_make_pam", "var2_make_pam"]
                            ].all(axis=1)
                            ind_cas_targ_pairs["both_break"] = ind_cas_targ_pairs[
                                ["var1_break_pam", "var2_break_pam"]
                            ].all(axis=1)
                            ind_cas_targ_pairs[
                                "one_make_one_break_1"
                            ] = ind_cas_targ_pairs[
                                ["var1_make_pam", "var2_break_pam"]
                            ].all(
                                axis=1
                            )
                            ind_cas_targ_pairs[
                                "one_make_one_break_2"
                            ] = ind_cas_targ_pairs[
                                ["var2_make_pam", "var1_break_pam"]
                            ].all(
                                axis=1
                            )
                            gens_replace = {
                                "0|1": "hap2",
                                "0|2": "hap2",
                                "0|3": "hap2",
                                "1|0": "hap1",
                                "2|0": "hap1",
                                "3|0": "hap1",
                                "0|0": "not_het",
                                "1|1": "not_het",
                            }
                            ind_gens = gens[["pos", ind]].replace(gens_replace).copy()
                            ind_cas_targ_pairs["var1_hap"] = pd.merge(
                                ind_cas_targ_pairs.copy(),
                                ind_gens.copy(),
                                left_on="var1",
                                right_on="pos",
                                how="left",
                            )[ind]
                            ind_cas_targ_pairs["var2_hap"] = pd.merge(
                                ind_cas_targ_pairs.copy(),
                                ind_gens.copy(),
                                left_on="var2",
                                right_on="pos",
                                how="left",
                            )[ind]
                            ind_cas_targ_pairs["same_hap"] = np.where(
                                ind_cas_targ_pairs["var1_hap"]
                                == ind_cas_targ_pairs["var2_hap"],
                                True,
                                False,
                            )
                            ind_cas_targ_pairs["not_same_hap"] = ~ind_cas_targ_pairs[
                                "same_hap"
                            ]
                            if (
                                ind_cas_targ_pairs[["both_make", "same_hap"]]
                                .all(axis=1)
                                .any()
                                or ind_cas_targ_pairs[["both_break", "same_hap"]]
                                .all(axis=1)
                                .any()
                            ):
                                ind_targ_cas.append(True)
                                continue
                            # check if pair where one makes, one breaks a PAM, and on different haplotypes
                            elif (
                                ind_cas_targ_pairs[
                                    ["one_make_one_break_1", "not_same_hap"]
                                ]
                                .all(axis=1)
                                .any()
                                or ind_cas_targ_pairs[
                                    ["one_make_one_break_2", "not_same_hap"]
                                ]
                                .all(axis=1)
                                .any()
                            ):
                                ind_targ_cas.append(True)
                                continue
                            # all possibilities exhausted, this person just isn't targetable at this gene
                            else:
                                ind_targ_cas.append(False)
                                continue

        finaltargcols.append(f"targ_{cas}")
        final_targ[f"targ_{cas}"] = ind_targ_cas

    # add column summarizing targetability across assessed Cas varieties

    final_targ["targ_all"] = final_targ[finaltargcols].any(axis=1)

    # HDF has issues with certain characters

    translated_gene_name = translate_gene_name(gene)

    if args["--exhaustive"]:
        exh_df = pd.concat(overall_exh_list, sort=False).drop_duplicates()
        exh_df.to_csv(f"{out_prefix}_exh.tsv", sep="\t", index=False)

    if args["-g"]:
        guides_out = pair_guides(args["--guides"], variant1, variant2)
        guides_out.to_csv(f"{out_prefix}pair_guides.tsv", sep="\t", index=False)

    # make list of genes that actually get written to HDF5
    with open(f"{out_prefix}genes_evaluated.txt", "a+") as f:
        f.write(f"{translated_gene_name}\n")
    # write gene dat to file
    final_targ.to_hdf(f"{out_prefix}.h5", "all", comp_level=9, complib="blosc")
    add_metadata(
        f"{out_prefix}.h5",
        args,
        os.path.basename(__file__),
        __version__,
        f"ExcisionFinder:{translated_gene_name}",
    )

    logging.info("Done!")


if __name__ == "__main__":
    arguments = docopt(__doc__, version=__version__)
    if arguments["-v"]:
        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s %(name)s:%(levelname)s ]%(message)s",
        )
    else:
        logging.basicConfig(
            level=logging.ERROR,
            format="[%(asctime)s %(name)s:%(levelname)s ]%(message)s",
        )
    logging.info(arguments)
    main(arguments)
