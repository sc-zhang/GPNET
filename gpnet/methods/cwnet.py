from gpnet.simulator.calc import calc_pheno
from gpnet.algorithm.SA import SA
from gpnet.io.data_io import DataLoader, DataSaver
from os import getpid
from numpy import sum, average, std, log


class CWNET:
    def __init__(self, in_file, weight_file, out_file):
        self.__in_file = in_file
        self.__weight_file = weight_file
        self.__out_file = out_file
        self.__inc_type_info = [0]

        self.__nodes = None
        self.__edges = None

    def generate_network(self, type_info, genotypes, phenotypes):
        converted_data = {}
        total_site_cnt = sum(type_info)

        for site in range(total_site_cnt):
            converted_data[site] = []
            for _ in range(len(genotypes)):
                if genotypes[_][site] == 1:
                    converted_data[site].append(phenotypes[_])

        self.__nodes = {}
        idx = 0
        for site in sorted(converted_data,
                           key=lambda x: [int(log(average(converted_data[x])))
                                          if len(converted_data[x]) >= 1 and average(converted_data[x]) != 0 else 0,
                                          -std(converted_data[x])
                                          if len(converted_data[x]) > 1 else 0]):
            self.__nodes[site] = idx
            idx += 1
        '''
        for site in converted_data:
            self.__nodes[site] = 0
            if len(converted_data[site]) >= 1:
                self.__nodes[site] += average(converted_data[site])
            if len(converted_data[site]) > 1:
                self.__nodes[site] -= 3 * std(converted_data[site])
        '''
        converted_data = {}

        for s1 in range(len(type_info)-1):
            for s2 in range(s1+1, len(type_info)):
                for t1 in range(type_info[s1]):
                    for t2 in range(type_info[s2]):
                        idx1 = self.__inc_type_info[s1]+t1
                        idx2 = self.__inc_type_info[s2]+t2

                        pair = tuple([idx1, idx2])
                        if pair not in converted_data:
                            converted_data[pair] = []
                        for _ in range(len(genotypes)):
                            genotype = genotypes[_]
                            if genotype[idx1] == genotype[idx2] == 1:
                                converted_data[pair].append(phenotypes[_])

        self.__edges = {}
        idx = 0
        for pair in sorted(converted_data,
                           key=lambda x: [int(log(average(converted_data[x])))
                                          if len(converted_data[x]) >= 1 and average(converted_data[x]) > 0 else 0,
                                          -std(converted_data[x])
                                          if len(converted_data[x]) > 1 else 0]):
            idx1, idx2 = pair
            if idx1 not in self.__edges:
                self.__edges[idx1] = {}
            self.__edges[idx1][idx2] = idx
            idx += 1
        '''
        for pair in converted_data:
            idx1, idx2 = pair
            if idx1 not in self.__edges:
                self.__edges[idx1] = {}
            self.__edges[idx1][idx2] = 0
            if len(converted_data[pair]) >= 1:
                self.__edges[idx1][idx2] += average(converted_data[pair])
            if len(converted_data[pair]) > 1:
                self.__edges[idx1][idx2] -= 3*std(converted_data[pair])
        '''

    def calc_score(self, data):
        score = 0
        available_sites = []
        for idx in range(len(data)):
            if data[idx] == 1:
                available_sites.append(idx)

        for idx in available_sites:
            score += self.__nodes[idx]

        for i in range(len(available_sites)-1):
            for j in range(i+1, len(available_sites)):
                score += self.__edges[available_sites[i]][available_sites[j]]
        return score

    def run(self):
        print("\tPID:%d Loading data" % getpid())
        dl = DataLoader()
        dl.load_genotype(self.__in_file)
        dl.load_weight_data(self.__weight_file)
        genotypes = dl.genotypes
        phenotypes = dl.phenotypes
        type_info = dl.type_info
        single_weight = dl.single_weight
        multi_weight = dl.multi_weight
        for type_cnt in type_info:
            self.__inc_type_info.append(self.__inc_type_info[-1] + type_cnt)

        print("\tPID:%d Generating network" % getpid())
        self.generate_network(type_info, genotypes, phenotypes)
        print("\tPID:%d Running SA" % getpid())
        sa = SA(type_info, self.calc_score, iterate=100, alpha=0.99)
        sa.run()
        best_sa_data = sa.data
        best_sa_pheno = calc_pheno(best_sa_data, single_weight, multi_weight)

        print("\tPID:%d Saving predict data" % getpid())
        ds = DataSaver(self.__out_file)
        ds.save_data(type_info, [best_sa_data], [best_sa_pheno])
