from numpy import random


def calc_pheno(genotype, single_weight, multi_weight):
    __RAND_EFFECT = [-100, 100]
    pheno = 0
    for _ in range(len(genotype)):
        pheno += genotype[_] * single_weight[_]
    external = 0
    matches = set()
    # high-order co-effect will mask low-order co-effect
    for _ in sorted(multi_weight, reverse=True):
        for sites in multi_weight[_]:
            is_match = True
            for site in sites:
                if genotype[site] == 0:
                    is_match = False
                    break
            if is_match:
                for already_matched in matches:
                    if len(sites) == len(set(sites).intersection(set(already_matched))):
                        is_match = False
                        break
            if is_match:
                external += multi_weight[_][sites]
                matches.add(sites)
    pheno += external

    # Simulate environment effect
    pheno += random.randint(__RAND_EFFECT[0], __RAND_EFFECT[1])
    return pheno
