from gpnet.simulator.calc import calc_pheno
from gpnet.algorithm.SA import SA
from gpnet.algorithm.GA import GA
from gpnet.io.data_io import DataLoader, DataSaver
from gpnet.io.message import Message
from os import getpid
from numpy import sum, average, std, isnan


class CWNET:
    def __init__(
        self,
        in_file,
        weight_file,
        is_lower_better,
        is_normalization,
        is_store_iter,
        top_edges,
        min_edge_cnt,
        step_size,
        optim_args,
        seed,
        out_file,
    ):
        self.__in_file = in_file
        self.__weight_file = weight_file
        self.__is_lower_better = is_lower_better
        self.__is_normalization = is_normalization
        self.__is_store_iter = is_store_iter
        self.__optim_args = optim_args
        self.__seed = seed
        self.__out_file = out_file
        self.__type_info = None
        self.__allele_name = None

        self.__nodes = None
        self.__edges = None

    def __convert_data(self, data_list):
        converted_data_db = {}
        if len(data_list) == 0:
            return converted_data_db

        min_avg = data_list[0][1]
        max_avg = data_list[0][1]

        min_std = data_list[0][2]
        max_std = data_list[0][2]

        for _, avg, stdv in data_list:
            if not isnan(avg):
                if avg < min_avg or isnan(min_avg):
                    min_avg = avg
                if avg > max_avg or isnan(max_avg):
                    max_avg = avg
            if not isnan(stdv):
                if stdv < min_std or isnan(min_std):
                    min_std = stdv
                if stdv > max_std or isnan(max_std):
                    max_std = stdv

        mm_avg = max_avg - min_avg
        mm_std = max_std - min_std
        if mm_avg == 0:
            mm_avg = 1.0
        if mm_std == 0:
            mm_std = 1.0
        for _ in range(len(data_list)):
            if self.__is_lower_better:
                if isnan(data_list[_][1]):
                    data_list[_][1] = max_avg
                if isnan(data_list[_][2]):
                    data_list[_][2] = max_std
            else:
                if isnan(data_list[_][1]):
                    data_list[_][1] = min_avg
                if isnan(data_list[_][2]):
                    data_list[_][2] = min_std
            if self.__is_normalization:
                data_list[_][1] = int((data_list[_][1] - min_avg) * 99.0 / mm_avg)
                data_list[_][2] = int((data_list[_][2] - min_std) * 99.0 / mm_std)

        # average higher and stdev lower is better
        for _, avg, stdv in sorted(data_list, key=lambda x: [x[1], -x[2]]):
            if isinstance(_, int):
                if self.__is_lower_better:
                    converted_data_db[_] = (
                        avg + stdv / 100.0 if self.__is_normalization else avg
                    )
                else:
                    converted_data_db[_] = (
                        avg + (0.99 - stdv / 100.0) if self.__is_normalization else avg
                    )
            else:
                idx1, idx2 = _
                if idx1 not in converted_data_db:
                    converted_data_db[idx1] = {}
                if self.__is_lower_better:
                    converted_data_db[idx1][idx2] = (
                        avg + stdv / 100.0 if self.__is_normalization else avg
                    )
                else:
                    converted_data_db[idx1][idx2] = (
                        avg + (0.99 - stdv / 100.0) if self.__is_normalization else avg
                    )
        return converted_data_db

    def generate_network(self, genotypes, phenotypes):
        converted_data = {}
        total_site_cnt = sum(self.__type_info)

        for site in range(total_site_cnt):
            converted_data[site] = []
            for _ in range(len(genotypes)):
                if genotypes[_][site] == 1:
                    converted_data[site].append(phenotypes[_])

        node_list = [
            [
                _,
                (
                    average(converted_data[_])
                    if len(converted_data[_]) >= 1
                    else float("nan")
                ),
                std(converted_data[_]) if len(converted_data[_]) >= 1 else float("nan"),
            ]
            for _ in range(total_site_cnt)
        ]

        self.__nodes = self.__convert_data(node_list)

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

        edge_list = [
            [
                pair,
                (
                    average(converted_data[pair])
                    if len(converted_data[pair]) > 0
                    else float("nan")
                ),
                (
                    std(converted_data[pair])
                    if len(converted_data[pair]) > 1
                    else float("nan")
                ),
            ]
            for pair in converted_data
        ]

        self.__edges = self.__convert_data(edge_list)

    def predict(self, data):
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
        """
        for idx1 in range(len(self.__allele_name) - 1):
            for idx2 in range(idx1 + 1, len(self.__allele_name)):
        """
        for idx1 in sorted(self.__edges):
            for idx2 in sorted(self.__edges[idx1]):
                edge_weight.append(
                    [
                        self.__allele_name[idx1],
                        self.__allele_name[idx2],
                        self.__edges[idx1][idx2],
                    ]
                )

        # generate additional information of weight of nodes and weight of edges, sorted by weight descend
        additional_info = ["#\n# Nodes weights"]
        for _ in sorted(
            node_weight, key=lambda x: x[-1], reverse=(not self.__is_lower_better)
        ):
            additional_info.append("# %s" % (" ".join(map(str, _))))

        additional_info.append("#\n# Edges weights")
        for _ in sorted(
            edge_weight, key=lambda x: x[-1], reverse=(not self.__is_lower_better)
        ):
            additional_info.append("# %s" % (" ".join(map(str, _))))

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

        Message.info("\tPID:%d Running Optimization" % getpid())
        if self.__optim_args["Method"] == "SA":
            sa = SA(
                self.__type_info,
                self.__allele_name,
                self.__out_file,
                self.predict,
                self.__is_lower_better,
                self.__is_store_iter,
                iterate=self.__optim_args["iter_cnt"],
                alpha=self.__optim_args["iter_alpha"],
                seed=self.__seed,
            )
            sa.run()
            best_data = sa.data
        elif self.__optim_args["Method"] == "GA":
            ga = GA(
                self.__type_info,
                self.__allele_name,
                self.__out_file,
                self.predict,
                self.__is_lower_better,
                self.__is_store_iter,
                n_generations=self.__optim_args["n_generations"],
                pop_size=self.__optim_args["pop_size"],
                cross_rate=self.__optim_args["cross_rate"],
                mutation_rate=self.__optim_args["mutation_rate"],
                seed=self.__seed,
            )
            ga.run()
            best_data = ga.data

        # best_pheno can only be calculated while data is generated by gpnet.py simulator
        if self.__weight_file:
            dl.load_weight_data(self.__weight_file)
            single_weight = dl.single_weight
            multi_weight = dl.multi_weight
            best_pheno = calc_pheno(best_data, single_weight, multi_weight)
        else:
            best_pheno = "Not support"

        Message.info("\tPID:%d Saving predict data" % getpid())
        ds = DataSaver(self.__out_file)
        ds.save_data(
            self.__type_info,
            self.__allele_name,
            [best_data],
            [best_pheno],
            self._get_node_edge_weight(),
        )
        Message.info("\tPID:%d Saved" % getpid())
