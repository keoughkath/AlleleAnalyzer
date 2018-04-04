# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script uses set cover theory to optimize sgRNA pair choice such that you get either the least number of probe pairs
to cover the most people. This can be done by specifying either the maximum number of probes you want to use or the
minimum number of people you would like to cover with the returned probe set. Kathleen Keough & Svetlana Lyalina 2017-2018.

Usage:
        optimize_ppl_covered.py --type=<type> <mp> <infile> <outprefix> [--guides] <guides>

Arguments:
        type                type of analysis, either max_probes or min_prop
        mp                  if type = max_probes, this is the maximum number of sgRNA pairs to consider,
                            if type = min_prop, this is the minimum proportion of the population that must be covered
        infile              this is the dataframe for the gene or locus being analyzed by gen_arcplot_df.py
        outprefix           where to save the outputted files, which include one file for the individuals covered
                            and one file for the sgRNA pairs used
Options:
        --guides            Output guides as well for identified variant pairs. 
        guides              Must enter this file, generate by gen_sgRNAs.py, if --guides is specified
"""
import numpy as np
import pandas as pd
from pulp import *
import os
from docopt import docopt


def optimize_probes(probes_to_people,
                    min_prop_covered=None,
                    max_probes=None):
    num_probes = probes_to_people.shape[0] 
    num_people = probes_to_people.shape[1]

    if min_prop_covered is None and max_probes is not None:
        problem_type = "maximize coverage"
    elif max_probes is None and min_prop_covered is not None:
        problem_type = "minimize probes"
    else:
        raise ValueError("Must provide either minimum fraction of coverage or maximum number of probes")

    if problem_type == "minimize probes":
        prob = LpProblem("Solving for minimal number of probes needed to cover " + str(100*min_prop_covered) + "% of people",LpMinimize)
    elif problem_type == "maximize coverage":
        prob = LpProblem("Solving for maximal cover with " + str(max_probes) + " max allowed",LpMaximize)
    probes = [LpVariable("g{}".format(i+1), cat="Binary") for i in range(num_probes)]
    people = [LpVariable("s{}".format(i+1), cat="Binary") for i in range(num_people)]

    # Add objective
    # A) minimize number of probes used
    if problem_type == "minimize probes":
        total_probes = sum(probes)
        prob += total_probes
    # B) maximize coverage
    elif problem_type == "maximize coverage":
        total_people = sum(people)
        prob += total_people

    # Add constraint
    # A) must cover at least x% of people
    if problem_type == "minimize probes":
        total_people = sum(people)
        prob += total_people >= (min_prop_covered*num_people)
    # B) must use at most x probes
    elif problem_type == "maximize coverage":
        total_probes = sum(probes)
        prob += total_probes <= max_probes

    # Add constraints arising from the probes to people mapping
    for idx, person in enumerate(people):
        #ugh that indexing is so ugly, should make better
        prob += sum([probes[i] for i in np.where(probes_to_people[:,idx])[0].tolist()]) - person >= 0

    if GUROBI().available():
        status = prob.solve(GUROBI())
    elif GLPK().available():
        status = prob.solve(GLPK())
    else:
        status = prob.solve()

    return({
        'probe usage status': [probe.value() for probe in probes],
        'person covered status': [person.value() for person in people],
        'solution status': status 
        })


def get_people(solution, samples):
    """
    Get people covered by solution.
    :param solution: outdict from optimize_probes
    :param samples: All sample available to be covered.
    :return: The individuals covered by the solution.
    """
    mask = pd.Series(solution['person covered status']).astype(bool)
    people_covered = pd.Series(samples)[mask].tolist()
    return (people_covered)


def get_pairs(solution, targ_pairs):
    """
    Get pairs used.
    :param solution: outdict from optimize_probes
    :param targ_pairs: Pandas dataframe of  pairs with the individuals they can target allele-specifically.
    :return: The pairs used.
    """
    mask = pd.Series(solution['probe usage status']).astype(bool)
    pairs_used = targ_pairs.loc[mask.astype(bool)][['var1', 'var2']]
    return (pairs_used)

def main(args):
    run_type = args['--type']
    mp = args['<mp>']
    indf = pd.read_csv(args['<infile>'], sep='\t')
    outprefix = args['<outprefix>']
    if run_type == 'max_probes':
        outdict = optimize_probes(indf.iloc[:,2:].applymap(int).as_matrix(), max_probes=int(mp))
    elif run_type == 'min_prop':
        outdict = optimize_probes(indf.iloc[:,2:].applymap(int).as_matrix(), min_prop_covered=float(mp))
    else:
        print('run_type must be either max_probes or min_prop.')
        exit()
    ppl_covered = get_people(outdict, indf.columns[2:])
    pairs_used = get_pairs(outdict, indf)
    with open(outprefix + '_ppl_covered.txt', 'w') as f:
        for person in ppl_covered:
            f.write(person + '\n')
    if args['--guides']:
        if not args['<guides>']:
            print('Must provides <guides> file. Not outputting guides.')
        else:
            guides_df = pd.read_csv(args['<guides>'], sep='\t')
            guide_pairs_out = pd.DataFrame()
            guide_pairs_out['variant1'] = pairs_used['var1']
            guide_pairs_out['variant2'] = pairs_used['var2']
            guide_pairs_out['variant_position'] = guide_pairs_out['variant1']
            guide_pairs_out = guide_pairs_out.merge(guides_df, how='left', on=['variant_position'])
            guide_pairs_out['gRNA_ref_guide_1'] = guide_pairs_out['gRNA_ref']
            guide_pairs_out['gRNA_alt_guide_1'] = guide_pairs_out['gRNA_alt']
            guide_pairs_out['ref_g1'] = guide_pairs_out['ref']
            guide_pairs_out['alt_g1'] = guide_pairs_out['alt']
            guide_pairs_out['variant_position_in_guide_g1'] = guide_pairs_out['variant_position_in_guide']
            guide_pairs_out = guide_pairs_out[['variant1','variant2','gRNA_ref_guide_1','gRNA_alt_guide_1',
            'ref_g1','alt_g1','variant_position_in_guide_g1']]
            guide_pairs_out['variant_position'] = guide_pairs_out['variant2']
            guide_pairs_out = guide_pairs_out.merge(guides_df, how='left', on=['variant_position'])
            guide_pairs_out['gRNA_ref_guide_2'] = guide_pairs_out['gRNA_ref']
            guide_pairs_out['gRNA_alt_guide_2'] = guide_pairs_out['gRNA_alt']
            guide_pairs_out['ref_g2'] = guide_pairs_out['ref']
            guide_pairs_out['alt_g2'] = guide_pairs_out['alt']
            guide_pairs_out['variant_position_in_guide_g2'] = guide_pairs_out['variant_position_in_guide']
            guides_out = guide_pairs_out[['variant1','variant2','gRNA_ref_guide_1','gRNA_alt_guide_1',
            'ref_g1','alt_g1','variant_position_in_guide_g1','gRNA_ref_guide_2','gRNA_alt_guide_2',
            'ref_g2','alt_g2','variant_position_in_guide_g2']].drop_duplicates()
            guides_out.to_csv(f'{outprefix}_optimized_pair_guides.tsv', sep='\t', index=False)
    pairs_used.to_csv(outprefix + '_pairs_used.txt', sep='\t', index=False)

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    main(arguments)
