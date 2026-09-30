from numpy import array, random, sum
from wpnet_module.gpnet.simulator.calc import calc_pheno


class Simulator:
    def __init__(
        self, site_cnt, avg_type_cnt, noise_ratio, sample_cnt, is_hybrid, pat_cnt
    ):
        random.seed()
        self.__site_cnt = site_cnt
        self.__avg_type_cnt = avg_type_cnt
        self.__noise_ratio = noise_ratio
        self.__sample_cnt = sample_cnt
        self.__is_hybrid = is_hybrid
        self.__pat_cnt = pat_cnt
        self.__total_type_cnt = 0

        self.__NOISE_EFFECT = [0, 20]
        self.__NORMAL_EFFECT = [21, 200]

        self.genotypes = None
        self.type_info = None
        self.single_weight = None
        self.multi_weight = None
        self.phenotypes = None

    def sim_genotypes(self):
        # random set variant type count for each variant site
        self.type_info = []

        """
        for _ in range(self.__site_cnt):
            self.type_info.append(random.randint(1, self.__avg_type_cnt * 2 + 2))
        """
        self.type_info = random.binomial(
            self.__avg_type_cnt * 2 + 1, 0.5, self.__site_cnt
        )

        max_type_cnt = self.__avg_type_cnt * 2 + 1
        for _ in range(len(self.type_info)):
            if self.type_info[_] > max_type_cnt:
                self.type_info[_] = max_type_cnt
            elif self.type_info[_] <= 0:
                self.type_info[_] = 1

        self.__total_type_cnt = sum(self.type_info)

        inc_list = [0]
        for type_cnt in self.type_info:
            inc_list.append(inc_list[-1] + type_cnt)

        """
        self.genotypes = array([[0 for __ in range(self.__total_type_cnt)] for _ in range(self.__sample_cnt)])
        for smp_idx in range(self.__sample_cnt):
            for var_idx in range(self.__site_cnt):
                self.genotypes[smp_idx][random.randint(0, self.type_info[var_idx]) + inc_list[var_idx]] = 1
                random_idx = random.poisson(lam=self.type_info[var_idx]/1.4)
                if random_idx < 0 or random_idx >= self.type_info[var_idx]:
                    random_idx = 0
                self.genotypes[smp_idx][random_idx + inc_list[var_idx]] = 1
        """

        self.genotypes = []
        if self.__is_hybrid:
            pat_pool = set()
            for pat_idx in range(self.__pat_cnt):
                pat_geno = [0 for _ in range(self.__total_type_cnt)]
                for var_idx in range(self.__site_cnt):
                    pat_geno[
                        random.randint(0, self.type_info[var_idx]) + inc_list[var_idx]
                    ] = 1

                # avoid duplicate samples, if the choise less than sample count, it may cause deadloop, for that
                # test_cnt is set as test times, max test should not over 10000 times.
                test_cnt = 0
                while tuple(pat_geno) in pat_pool and test_cnt < 10000:
                    pat_geno = [0 for _ in range(self.__total_type_cnt)]
                    for var_idx in range(self.__site_cnt):
                        pat_geno[
                            random.randint(0, self.type_info[var_idx])
                            + inc_list[var_idx]
                        ] = 1
                    test_cnt += 1

                pat_pool.add(tuple(pat_geno))
            pat_pool = list(pat_pool)
            for smp_idx in range(self.__sample_cnt):
                random.shuffle(pat_pool)
                genotype = [0 for _ in range(self.__total_type_cnt)]
                pat_geno = pat_pool[0]
                mat_geno = pat_pool[1]
                for var_idx in range(self.__site_cnt):
                    if random.rand() < 0.5:
                        for _ in range(inc_list[var_idx], inc_list[var_idx + 1]):
                            genotype[_] = pat_geno[_]
                    else:
                        for _ in range(inc_list[var_idx], inc_list[var_idx + 1]):
                            genotype[_] = mat_geno[_]
                self.genotypes.append(genotype)
        else:
            genotype_set = set()
            for smp_idx in range(self.__sample_cnt):
                genotype = [0 for _ in range(self.__total_type_cnt)]
                for var_idx in range(self.__site_cnt):
                    genotype[
                        random.randint(0, self.type_info[var_idx]) + inc_list[var_idx]
                    ] = 1

                # avoid duplicate samples, if the choise less than sample count, it may cause deadloop, for that
                # test_cnt is set as test times, max test should not over 10000 times.
                test_cnt = 0
                while tuple(genotype) in genotype_set and test_cnt < 10000:
                    genotype = [0 for _ in range(self.__total_type_cnt)]
                    for var_idx in range(self.__site_cnt):
                        genotype[
                            random.randint(0, self.type_info[var_idx])
                            + inc_list[var_idx]
                        ] = 1
                    test_cnt += 1

                genotype_set.add(tuple(genotype))
                self.genotypes.append(genotype)
        self.genotypes = array(self.genotypes)

    def sim_site_effects(self):
        # simulate single site effect
        self.single_weight = array([0 for _ in range(self.__total_type_cnt)])

        # for each site, if variant type more than 1, the last variant type is set as lost,
        # the effect of 'lost' variant is also set to 0
        noise_sites = set()
        site_effect_idx = 0
        for var_idx in range(self.__site_cnt):
            for type_idx in range(self.type_info[var_idx]):
                if (
                    type_idx == self.type_info[var_idx] - 1
                    and self.type_info[var_idx] > 1
                ):
                    self.single_weight[site_effect_idx] = 0
                else:
                    sign = 1 if random.rand() < 0.5 else -1
                    if random.rand() < self.__noise_ratio:
                        self.single_weight[site_effect_idx] = sign * random.randint(
                            self.__NOISE_EFFECT[0], self.__NOISE_EFFECT[1]
                        )
                        """
                        self.single_weight[site_effect_idx] = random.normal(loc=0, scale=self.__NOISE_EFFECT[1])
                        """
                        noise_sites.add(site_effect_idx)
                    else:

                        self.single_weight[site_effect_idx] = sign * random.randint(
                            self.__NORMAL_EFFECT[0], self.__NORMAL_EFFECT[1]
                        )
                        """
                        self.single_weight[site_effect_idx] = random.normal(loc=0, scale=self.__NORMAL_EFFECT[1])
                        """

                site_effect_idx += 1

        # simulate multi sites co-effect
        self.multi_weight = {}
        # simulate 2 to 9 sites co-effect
        for co_site_cnt in range(2, 10):
            self.multi_weight[co_site_cnt] = {}
            lower = int(self.__total_type_cnt * 0.1)
            upper = int(self.__total_type_cnt * 0.25 + 1)
            comb_cnt = random.randint(lower, upper)

            while len(self.multi_weight[co_site_cnt]) < comb_cnt:
                # co-effect only appear among different sites
                tmp_sites = set()
                while len(tmp_sites) < co_site_cnt:
                    tmp_sites.add(random.randint(0, self.__site_cnt))

                # for each site, if the variant type lager than 1,
                # the last type is set as lost, means it must not appear in co-effect sites
                tmp_site_with_types = []

                for site in tmp_sites:
                    # the co-effect must not appear while one site is marked as lost
                    site_with_type_idx = site + random.randint(
                        0, 1 if self.type_info[site] == 1 else self.type_info[site] - 1
                    )

                    tmp_site_with_types.append(site_with_type_idx)

                # simulate additive effect and epistatic effect
                if random.random() < 1.0 / (100 * (10 - co_site_cnt)):
                    effect_type = "epi"
                else:
                    effect_type = "add"

                curr_effect = 0
                if effect_type == "add":
                    for _ in tmp_site_with_types:
                        if _ not in noise_sites:
                            curr_effect += self.single_weight[_]
                else:
                    for _ in tmp_site_with_types:
                        if _ not in noise_sites:
                            if curr_effect < self.single_weight[_]:
                                curr_effect = self.single_weight[_]
                self.multi_weight[co_site_cnt][
                    tuple(sorted(tmp_site_with_types))
                ] = curr_effect

    def sim_phenotypes(self):
        self.phenotypes = array([0 for _ in range(self.__sample_cnt)])
        idx = 0
        for genotype in self.genotypes:
            self.phenotypes[idx] = calc_pheno(
                genotype, self.single_weight, self.multi_weight
            )
            idx += 1
