from gpnet.simulator.calc import calc_pheno
from gpnet.algorithm.SA import SA
from gpnet.io.data_io import DataLoader, DataSaver
from gpnet.io.message import Message
from os import getpid
from numpy import sum, average, std


class CWNET:
    def __init__(self, in_file, weight_file, is_lower_better, out_file):
        self.__in_file = in_file
        self.__weight_file = weight_file
        self.__is_lower_better = is_lower_better
        self.__out_file = out_file
        self.__type_info = None
        self.__allele_name = None

        self.__nodes = None
        self.__edges = None

    def generate_network(self, genotypes, phenotypes):
        converted_data = {}
        total_site_cnt = sum(self.__type_info)

        for site in range(total_site_cnt):
            converted_data[site] = []
            for _ in range(len(genotypes)):
                if genotypes[_][site] == 1:
                    converted_data[site].append(phenotypes[_])
        if self.__is_lower_better:
            unordered_list = [[_,
                               average(converted_data[_]) if len(converted_data[_]) >= 1 else 100,
                               std(converted_data[_]) if len(converted_data[_]) >= 1 else 100]
                              for _ in range(total_site_cnt)]
        else:
            unordered_list = [[_,
                               average(converted_data[_]) if len(converted_data[_]) >= 1 else 0,
                               std(converted_data[_]) if len(converted_data[_]) >= 1 else 0]
                              for _ in range(total_site_cnt)]
        '''
        self.__nodes = {}

        avg_list = [average(converted_data[_]) if len(converted_data[_]) >= 1 else 0 for _ in range(total_site_cnt)]
        min_avg = min(avg_list)
        max_avg = max(avg_list)
        # avoid divide zero
        mm = max_avg - min_avg + 1
        avg_list = [int((_-min_avg)*100./mm) for _ in avg_list]
        std_list = [std(converted_data[_]) if len(converted_data[_]) >= 1 else 0 for _ in range(total_site_cnt)]
        min_std = min(std_list)
        max_std = max(std_list)
        mm = max_std - min_std + 1
        std_list = [int((_-min_avg)*100./mm) for _ in std_list]

        idx = 0
        for site in sorted(converted_data, key=lambda x: [avg_list[x], -std_list[x]]):
            self.__nodes[site] = idx
            idx += 1
        '''

        converted_data = {}
        inc_type_info = [0]
        for _ in self.__type_info:
            inc_type_info.append(inc_type_info[-1] + _)

        for site1 in range(len(self.__type_info) - 1):
            for site2 in range(site1 + 1, len(self.__type_info)):
                for type1 in range(self.__type_info[site1]):
                    for type2 in range(self.__type_info[site2]):
                        idx1 = inc_type_info[site1] + type1
                        idx2 = inc_type_info[site2] + type2
                        pair = tuple([idx1, idx2])
                        if pair not in converted_data:
                            converted_data[pair] = []
                        for _ in range(len(genotypes)):
                            genotype = genotypes[_]
                            if genotype[idx1] == genotype[idx2] == 1:
                                converted_data[pair].append(phenotypes[_])
        '''
        for idx1 in range(len(self.__allele_name) - 1):
            for idx2 in range(idx1 + 1, len(self.__allele_name)):
                pair = tuple([idx1, idx2])
                if pair not in converted_data:
                    converted_data[pair] = []
                for _ in range(len(genotypes)):
                    genotype = genotypes[_]
                    if genotype[idx1] == genotype[idx2] == 1:
                        converted_data[pair].append(phenotypes[_])
        '''

        for pair in converted_data:
            if self.__is_lower_better:
                unordered_list.append([pair,
                                       average(converted_data[pair] if len(converted_data[pair]) > 0 else 100),
                                       std(converted_data[pair]) if len(converted_data[pair]) > 1 else 100])
            else:
                unordered_list.append([pair,
                                       average(converted_data[pair] if len(converted_data[pair]) > 0 else 0),
                                       std(converted_data[pair]) if len(converted_data[pair]) > 1 else 0])
        '''
        self.__edges = {}

        avg_db = {pair: average(converted_data[pair]) if len(converted_data[pair]) > 0 else 0
                  for pair in converted_data}
        min_avg = min([avg_db[_] for _ in avg_db])
        max_avg = max([avg_db[_] for _ in avg_db])
        mm = max_avg - min_avg + 1
        avg_db = {pair: int((avg_db[pair]-min_avg)*100./mm) for pair in avg_db}
        
        std_db = {pair: std(converted_data[pair]) if len(converted_data[pair]) > 1 else 0
                  for pair in converted_data}
        min_std = min([std_db[_] for _ in std_db])
        max_std = max([std_db[_] for _ in std_db])
        mm = max_std - min_std + 1
        std_db = {pair: int((std_db[pair]-min_std)*100./mm) for pair in std_db}

        idx = 0
        for pair in sorted(converted_data, key=lambda x: [avg_db[x], -std_db[x]]):
            idx1, idx2 = pair
            if idx1 not in self.__edges:
                self.__edges[idx1] = {}
            self.__edges[idx1][idx2] = idx
            idx += 1
        '''

        min_avg = unordered_list[0][1]
        max_avg = unordered_list[0][1]

        min_std = unordered_list[0][2]
        max_std = unordered_list[0][2]

        for _, avg, stdv in unordered_list:
            if avg < min_avg:
                min_avg = avg
            if avg > max_avg:
                max_avg = avg
            if stdv < min_std:
                min_std = stdv
            if stdv > max_std:
                max_std = stdv

        mm_avg = max_avg - min_avg + 1
        mm_std = max_std - min_std + 1

        for _ in range(len(unordered_list)):
            unordered_list[_][1] = int((unordered_list[_][1] - min_avg) * 100. / mm_avg)
            unordered_list[_][2] = int((unordered_list[_][2] - min_std) * 100. / mm_std)

        self.__nodes = {}
        self.__edges = {}
        # average higher and stdev lower is better
        for _, avg, stdv in sorted(unordered_list, key=lambda x: [x[1], -x[2]]):
            if isinstance(_, int):
                if self.__is_lower_better:
                    self.__nodes[_] = avg + stdv / 100.
                else:
                    self.__nodes[_] = avg + (1 - stdv / 100.)
            else:
                idx1, idx2 = _
                if idx1 not in self.__edges:
                    self.__edges[idx1] = {}
                if self.__is_lower_better:
                    self.__edges[idx1][idx2] = avg + stdv / 100.
                else:
                    self.__edges[idx1][idx2] = avg + (1 - stdv / 100.)

    def calc_score(self, data):
        score = 0
        available_sites = []
        for idx in range(len(data)):
            if data[idx] == 1:
                available_sites.append(idx)

        for idx in available_sites:
            score += self.__nodes[idx]

        for i in range(len(available_sites) - 1):
            for j in range(i + 1, len(available_sites)):
                score += self.__edges[available_sites[i]][available_sites[j]]
        return score

    def _get_node_edge_weight(self):
        node_weight = []
        edge_weight = []

        for idx in range(len(self.__allele_name)):
            node_weight.append([self.__allele_name[idx], self.__nodes[idx]])
        '''
        for idx1 in range(len(self.__allele_name) - 1):
            for idx2 in range(idx1 + 1, len(self.__allele_name)):
        '''
        for idx1 in sorted(self.__edges):
            for idx2 in sorted(self.__edges[idx1]):
                edge_weight.append([self.__allele_name[idx1], self.__allele_name[idx2], self.__edges[idx1][idx2]])

        # generate additional information of weight of nodes and weight of edges, sorted by weight descend
        additional_info = ["#\n# Nodes weights"]
        for _ in sorted(node_weight, key=lambda x: x[-1], reverse=(not self.__is_lower_better)):
            additional_info.append("# %s" % (' '.join(map(str, _))))

        additional_info.append("#\n# Edges weights")
        for _ in sorted(edge_weight, key=lambda x: x[-1], reverse=(not self.__is_lower_better)):
            additional_info.append("# %s" % (' '.join(map(str, _))))

        return additional_info

    def run(self):
        Message.info("\tPID:%d Loading data" % getpid())
        dl = DataLoader()
        dl.load_genotype(self.__in_file)
        genotypes = dl.genotypes
        phenotypes = dl.phenotypes
        self.__type_info = dl.type_info
        self.__allele_name = dl.allele_name

        Message.info("\tPID:%d Generating network" % getpid())
        self.generate_network(genotypes, phenotypes)
        Message.info("\tPID:%d Running SA" % getpid())
        sa = SA(self.__type_info, self.__allele_name, self.__out_file,
                self.calc_score, self.__is_lower_better, iterate=100, alpha=0.99)
        sa.run()
        best_sa_data = sa.data

        # best_sa_pheno can only be calculated while data is generated by gpnet.py simulator
        if self.__weight_file:
            dl.load_weight_data(self.__weight_file)
            single_weight = dl.single_weight
            multi_weight = dl.multi_weight
            best_sa_pheno = calc_pheno(best_sa_data, single_weight, multi_weight)
        else:
            best_sa_pheno = "Not support"

        Message.info("\tPID:%d Saving predict data" % getpid())
        ds = DataSaver(self.__out_file)
        ds.save_data(self.__type_info, self.__allele_name, [best_sa_data], [best_sa_pheno],
                     self._get_node_edge_weight())
