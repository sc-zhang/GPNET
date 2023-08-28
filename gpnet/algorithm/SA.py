from copy import deepcopy
from numpy import random, exp, sum
from gpnet.io.data_io import DataSaver


class SA:
    def __init__(self, type_info, allele_name, outfile,
                 func, is_lower_better, iterate=100, t0=100, t_final=0.01, alpha=0.99):
        self.__type_info = type_info
        self.__allele_name = allele_name
        self.__outfile = outfile+".iter"
        self.__site_cnt = len(type_info)
        self.__inc_type_info = [0]
        for type_cnt in self.__type_info:
            self.__inc_type_info.append(self.__inc_type_info[-1]+type_cnt)
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
        change_site = random.randint(self.__site_cnt)
        change_type = random.randint(self.__type_info[change_site])

        for cur_type in range(self.__inc_type_info[change_site], self.__inc_type_info[change_site+1]):
            if cur_type-self.__inc_type_info[change_site] != change_type:
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
        while self.__t > self.__t_final:
            for _ in range(self.__iterate):
                f = self.__func(self.data)
                new_data = self.__generate_new()
                f_new = self.__func(new_data)
                if self.__metrospolis(f, f_new):
                    self.data = deepcopy(new_data)
                    if not self.__iter_data or self.data != self.__iter_data[-1]:
                        self.__iter_data.append(self.data)
                        self.__iter_score.append(f_new)

            self.__t = self.__t * self.__alpha

        ds = DataSaver(self.__outfile)
        ds.save_iter(self.__type_info, self.__allele_name, self.__iter_data, self.__iter_score)
