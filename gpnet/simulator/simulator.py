from numpy import array, random, sum
from gpnet.simulator.calc import calc_pheno


class Simulator:
    def __init__(self, site_cnt, max_type_cnt, noise_ratio, sample_cnt):
        random.seed()
        self.__site_cnt = site_cnt
        self.__max_type_cnt = max_type_cnt
        self.__noise_ratio = noise_ratio
        self.__sample_cnt = sample_cnt
        self.__total_type_cnt = 0

        self.__NOISE_EFFECT = [0, 20]
        self.__NORMAL_EFFECT = [50, 100]

        self.genotypes = None
        self.type_info = None
        self.single_weight = None
        self.multi_weight = None
        self.phenotypes = None

    def sim_genotypes(self):
        # random set variant type count for each variant site
        self.type_info = []
        '''
        for _ in range(self.__site_cnt):
            self.type_info.append(random.randint(1, self.__max_type_cnt + 2))
        '''
        self.type_info = random.poisson(lam=self.__max_type_cnt / 3., size=self.__site_cnt)

        for _ in range(len(self.type_info)):
            if self.type_info[_] > self.__max_type_cnt:
                self.type_info[_] = self.__max_type_cnt
            elif self.type_info[_] <= 0:
                self.type_info[_] = 1
        self.__total_type_cnt = sum(self.type_info)

        inc_list = [0]
        for type_cnt in self.type_info:
            inc_list.append(inc_list[-1] + type_cnt)

        self.genotypes = array([[0 for __ in range(self.__total_type_cnt)] for _ in range(self.__sample_cnt)])
        for smp_idx in range(self.__sample_cnt):
            for var_idx in range(self.__site_cnt):
                random_idx = random.poisson(lam=self.type_info[var_idx]/1.4)
                if random_idx < 0 or random_idx >= self.type_info[var_idx]:
                    random_idx = 0
                self.genotypes[smp_idx][random_idx + inc_list[var_idx]] = 1

    def sim_site_effects(self):
        # simulate single site effect
        self.single_weight = array([0 for _ in range(self.__total_type_cnt)])

        # for each site, if variant type more than 1, the last variant type is set as lost,
        # the effect of 'lost' variant is also set to 0
        noise_sites = set()
        site_effect_idx = 0
        for var_idx in range(self.__site_cnt):
            for type_idx in range(self.type_info[var_idx]):
                if type_idx == self.type_info[var_idx] - 1 and self.type_info[var_idx] > 1:
                    self.single_weight[site_effect_idx] = 0
                else:
                    if random.rand() < self.__noise_ratio:
                        self.single_weight[site_effect_idx] = random.randint(self.__NOISE_EFFECT[0],
                                                                             self.__NOISE_EFFECT[1])
                        noise_sites.add(site_effect_idx)
                    else:
                        self.single_weight[site_effect_idx] = random.randint(self.__NORMAL_EFFECT[0],
                                                                             self.__NORMAL_EFFECT[1])
                site_effect_idx += 1

        # simulate multi sites co-effect
        self.multi_weight = {}
        # simulate 2 to 9 sites co-effect
        for co_site_cnt in range(2, 10):
            self.multi_weight[co_site_cnt] = {}
            lower = self.__total_type_cnt * .1
            upper = self.__total_type_cnt * .25 + 1
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
                    site_with_type_idx = site + random.randint(0, 1 if self.type_info[site] == 1 else
                                                               self.type_info[site] - 1)

                    tmp_site_with_types.append(site_with_type_idx)

                # simulate additive effect and epistatic effect
                if random.random() < 1./co_site_cnt:
                    effect_type = "add"
                else:
                    effect_type = "epi"

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
                self.multi_weight[co_site_cnt][tuple(sorted(tmp_site_with_types))] = curr_effect

    def sim_phenotypes(self):
        self.phenotypes = array([0 for _ in range(self.__sample_cnt)])
        idx = 0
        for genotype in self.genotypes:
            self.phenotypes[idx] = calc_pheno(genotype, self.single_weight, self.multi_weight)
            idx += 1
