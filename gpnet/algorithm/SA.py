from copy import deepcopy
from numpy import random, exp, sum
from gpnet.io.data_io import DataSaver


class SA:
    def __init__(self, type_info, allele_name, outfile,
                 func, is_lower_better, iterate=100, t0=100, t_final=0.01, alpha=0.99):
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
        self.__iterate = iterate
        self.__t0 = t0
        self.__t_final = t_final
        self.__t = t0
        self.__alpha = alpha
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

            for cur_type in range(self.__inc_type_info[change_site], self.__inc_type_info[change_site + 1]):
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
        f = self.__func(self.data)
        while self.__t > self.__t_final:
            for _ in range(self.__iterate):
                new_data = self.__generate_new()
                f_new = self.__func(new_data)
                if self.__metrospolis(f, f_new):
                    f = f_new
                    self.data = deepcopy(new_data)
                    if not self.__iter_data or self.data != self.__iter_data[-1]:
                        self.__iter_data.append(self.data)
                        self.__iter_score.append(f_new)

            self.__t = self.__t * self.__alpha

        ds = DataSaver(self.__outfile)
        ds.save_iter(self.__type_info, self.__allele_name, self.__iter_data, self.__iter_score)


class SelectSA:
    def __init__(self, select_count, type_info, best_genotype, func,
                 is_lower_better, iterate=100, t0=100, t_final=0.01, alpha=0.99):
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

        self.data = [0 for _ in range(self.__site_cnt)]
        # set first select count of genes to 1
        self.__data_init__cnt = 0
        for _ in range(self.__site_cnt):
            if _ in self.__non_absence_site:
                self.data[_] = 1
                self.__data_init__cnt += 1
            if self.__data_init__cnt >= self.__select_count:
                break
        self.geno = [0 for _ in range(self.__total_type_cnt)]
        for _ in range(self.__select_count):
            self.geno[self.geno_db[_]] = 1

    def __generate_new(self):
        new_data = deepcopy(self.data)
        src_set = set()
        for _ in range(len(self.data)):
            if self.data[_] == 1:
                src_set.add(_)

        src_list = list(src_set)
        non_absence_site = list(self.__non_absence_site)
        for _ in range(1, random.randint(2, 10)):
            if len(src_set) >= len(non_absence_site):
                continue
            src_set = set(src_list)
            idx = random.randint(len(src_list))
            src_site = src_list[idx]
            tgt_idx = random.randint(len(non_absence_site))
            tgt_site = non_absence_site[tgt_idx]
            while tgt_site in src_set:
                tgt_idx = random.randint(len(non_absence_site))
                tgt_site = non_absence_site[tgt_idx]
            new_data[src_site] = 0
            new_data[tgt_site] = 1
            src_list.remove(src_site)
            src_list.append(tgt_site)
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
        if self.__data_init__cnt < self.__select_count:
            return False
        f = self.__func(self.geno)
        while self.__t > self.__t_final:
            for _ in range(self.__iterate):
                new_data = self.__generate_new()
                new_geno = [0 for _ in range(self.__total_type_cnt)]
                for _ in range(len(new_data)):
                    if new_data[_] == 1:
                        new_geno[self.geno_db[_]] = 1
                f_new = self.__func(new_geno)
                if self.__metrospolis(f, f_new):
                    f = f_new
                    self.data = deepcopy(new_data)
                    self.geno = deepcopy(new_geno)

            self.__t = self.__t * self.__alpha
        return True
