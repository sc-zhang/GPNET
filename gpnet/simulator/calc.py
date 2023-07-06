from numpy import random


def calc_pheno(genotype, single_weight, multi_weight):
    __RAND_EFFECT = [0, 100]

    multi_site_effect = 0
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
                multi_site_effect += multi_weight[_][sites]
                matches.add(sites)

    # single site effect will be masked by co-effect
    single_site_effect = 0
    for _ in range(len(genotype)):
        if _ not in matches:
            single_site_effect += genotype[_] * single_weight[_]

    # Simulate environment effect
    pheno = single_site_effect + multi_site_effect + random.randint(__RAND_EFFECT[0], __RAND_EFFECT[1])
    return pheno
