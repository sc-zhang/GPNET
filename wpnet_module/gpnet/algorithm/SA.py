from copy import deepcopy
from numpy import random, exp, sum
from gpnet.io.data_io import DataSaver, DataLoader
from gpnet.algorithm import check_comb_exists


class SA:
    def __init__(
        self,
        type_info,
        allele_name,
        outfile,
        func,
        is_lower_better,
        is_store_iter,
        seed=None,
        **kwargs,
    ):
        """
        iterate: int, number of iterations, default=100
        t0: float, initial temperature, default=100.
        t_final: float, final temperature, default=0.01
        alpha: float, temperature decend rate, default=0.99
        """
        if seed:
            random.seed(seed)
        else:
            random.seed()
        self.__type_info = type_info
        self.__allele_name = allele_name
        self.__outfile = outfile + ".iter"
        self.__site_cnt = len(type_info)
        self.__inc_type_info = [0]
        for type_cnt in self.__type_info:
            self.__inc_type_info.append(self.__inc_type_info[-1] + type_cnt)
        self.__func = func
        self.__is_lower_better = is_lower_better
        self.__is_store_iter = is_store_iter
        self.__iterate = kwargs.get("iterate", 100)
        self.__t0 = kwargs.get("t0", 100.0)
        self.__t_final = kwargs.get("t_final", 0.01)
        self.__alpha = kwargs.get("alpha", 0.99)
        self.__t = self.__t0
        self.__total_type_cnt = sum(type_info)
        self.data = [0 for _ in range(self.__total_type_cnt)]
        self.__iter_data = []
        self.__iter_score = []
        # default type is first one
        for _ in self.__inc_type_info[:-1]:
            self.data[_] = 1

    def __generate_new(self):
        new_data = deepcopy(self.data)
        # random change most 10 site
        for _ in range(1, random.randint(2, 10)):
            change_site = random.randint(self.__site_cnt)
            change_type = random.randint(self.__type_info[change_site])

            for cur_type in range(
                self.__inc_type_info[change_site], self.__inc_type_info[change_site + 1]
            ):
                if cur_type - self.__inc_type_info[change_site] != change_type:
                    new_data[cur_type] = 0
                else:
                    new_data[cur_type] = 1

        return new_data

    def __metrospolis(self, f, f_new):
        if self.__is_lower_better:
            if f_new <= f:
                return 1
            else:
                p = exp((f - f_new) / self.__t)
            if random.rand() < p:
                return 1
            else:
                return 0
        else:
            if f_new >= f:
                return 1
            else:
                p = exp((f_new - f) / self.__t)
                if random.rand() < p:
                    return 1
                else:
                    return 0

    def run(self):
        best_data = deepcopy(self.data)
        f = self.__func(self.data)
        best_f = f
        while self.__t > self.__t_final:
            for _ in range(self.__iterate):
                new_data = self.__generate_new()
                f_new = self.__func(new_data)
                if self.__metrospolis(f, f_new):
                    f = f_new
                    self.data = deepcopy(new_data)
                    if self.__is_store_iter and (
                        not self.__iter_data or self.data != self.__iter_data[-1]
                    ):
                        self.__iter_data.append(self.data)
                        self.__iter_score.append(f_new)
                if f_new > best_f:
                    best_f = f_new
                    best_data = deepcopy(new_data)
            self.__t = self.__t * self.__alpha
        self.data = deepcopy(best_data)
        if self.__is_store_iter:
            ds = DataSaver(self.__outfile)
            ds.save_iter(
                self.__type_info,
                self.__allele_name,
                self.__iter_data,
                self.__iter_score,
            )


class SelectSA:
    def __init__(
        self,
        select_count,
        type_info,
        best_genotype,
        mat_file,
        func,
        is_lower_better,
        iterate=100,
        t0=100,
        t_final=0.01,
        alpha=0.99,
        seed=None,
    ):
        if seed:
            random.seed(seed)
        else:
            random.seed()
        self.__type_info = type_info
        self.__site_cnt = len(type_info)
        self.__inc_type_info = [0]
        for type_cnt in self.__type_info:
            self.__inc_type_info.append(self.__inc_type_info[-1] + type_cnt)
        self.__non_absence_site = set()
        self.__func = func
        self.__is_lower_better = is_lower_better
        self.__iterate = iterate
        self.__t0 = t0
        self.__t_final = t_final
        self.__t = t0
        self.__alpha = alpha
        self.__total_type_cnt = sum(type_info)
        self.__select_count = select_count
        self.geno_db = {}
        for _ in range(self.__site_cnt):
            sp = self.__inc_type_info[_]
            ep = self.__inc_type_info[_ + 1]
            for __ in range(sp, ep):
                if best_genotype[__] == 1:
                    if __ != ep - 1:
                        self.__non_absence_site.add(_)
                    self.geno_db[_] = __
                    break

        self.__col_masks = None
        self.__col_popcnt = None
        init_data = []
        if mat_file is not None:
            dl = DataLoader()
            dl.load_genotype(mat_file)
            pop_geno = dl.genotypes
            self.__col_masks, self.__col_popcnt = check_comb_exists.build_column_masks(
                pop_geno
            )
            for smp_idx in range(len(dl.genotypes)):
                init_data = []
                for _ in range(self.__site_cnt):
                    sp = self.__inc_type_info[_]
                    ep = self.__inc_type_info[_ + 1]
                    for __ in range(sp, ep):
                        if (
                            __ not in self.__non_absence_site
                            and best_genotype[__] == dl.genotypes[smp_idx][__] == 1
                        ):
                            init_data.append(_)
                            break
                    if len(init_data) >= self.__select_count:
                        break
                if len(init_data) >= self.__select_count:
                    break

        self.data = [0 for _ in range(self.__site_cnt)]
        # set first select count of genes to 1
        self.__data_init_cnt = 0
        if init_data:
            for _ in init_data:
                self.data[_] = 1
                self.__data_init_cnt += 1
                if self.__data_init_cnt >= self.__select_count:
                    break
        else:
            for _ in range(self.__site_cnt):
                if _ in self.__non_absence_site:
                    self.data[_] = 1
                    self.__data_init_cnt += 1
                if self.__data_init_cnt >= self.__select_count:
                    break
        self.geno = [0 for _ in range(self.__total_type_cnt)]
        for _ in range(self.__select_count):
            self.geno[self.geno_db[_]] = 1

    def __to_geno(self, new_data):
        new_geno = [0 for _ in range(self.__total_type_cnt)]
        new_idx = []
        for _ in range(len(new_data)):
            if new_data[_] == 1:
                new_geno[self.geno_db[_]] = 1
                new_idx.append(self.geno_db[_])
        return new_idx, new_geno

    def __generate_new(self):
        new_data = deepcopy(self.data)
        src_set = set()
        for _ in range(len(self.data)):
            if self.data[_] == 1:
                src_set.add(_)

        src_list = list(src_set)
        non_absence_site = list(self.__non_absence_site)

        idx = random.randint(len(src_list))
        src_site = src_list[idx]
        try_cnt = 0
        while try_cnt < 1e3:
            tgt_idx = random.randint(len(non_absence_site))
            tgt_site = non_absence_site[tgt_idx]
            while tgt_site in src_set:
                tgt_idx = random.randint(len(non_absence_site))
                tgt_site = non_absence_site[tgt_idx]
            new_data[src_site] = 0
            new_data[tgt_site] = 1
            new_idx, new_geno = self.__to_geno(new_data)
            if check_comb_exists.exists_covering_row(
                new_idx, self.__col_masks, self.__col_popcnt
            ):
                break
            # restore change
            new_data[src_site] = 1
            new_data[tgt_site] = 0
            try_cnt += 1
        return new_data, new_geno

    def __metrospolis(self, f, f_new):
        if self.__is_lower_better:
            if f_new <= f:
                return 1
            else:
                p = exp((f - f_new) / self.__t)
            if random.rand() < p:
                return 1
            else:
                return 0
        else:
            if f_new >= f:
                return 1
            else:
                p = exp((f_new - f) / self.__t)
                if random.rand() < p:
                    return 1
                else:
                    return 0

    def run(self):
        if self.__data_init_cnt < self.__select_count:
            return False
        best_geno = deepcopy(self.geno)
        best_data = deepcopy(self.data)
        f = self.__func(self.geno)
        best_f = f
        while self.__t > self.__t_final:
            for _ in range(self.__iterate):
                new_data, new_geno = self.__generate_new()
                f_new = self.__func(new_geno)
                if self.__metrospolis(f, f_new):
                    f = f_new
                    self.data = deepcopy(new_data)
                    self.geno = deepcopy(new_geno)
                if f_new > best_f:
                    best_f = f_new
                    best_data = deepcopy(new_data)
                    best_geno = deepcopy(new_geno)
            self.__t = self.__t * self.__alpha
        self.data = deepcopy(best_data)
        self.geno = deepcopy(best_geno)
        return True
