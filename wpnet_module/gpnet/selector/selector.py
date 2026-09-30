from wpnet_module.gpnet.algorithm.SA import SelectSA
import numpy as np
from wpnet_module.gpnet import fast_ridge


class SelectNet:
    def __init__(
        self,
        select_count,
        type_info,
        allele_name,
        model_file,
        mat_file,
        outfile,
        best_genotype,
        nodes,
        edges,
        is_lower_better,
        iter,
        alpha,
        seed,
    ):
        self.__type_info = type_info
        self.__allele_name = allele_name
        if model_file is not None:
            self.__model = fast_ridge.GraphModel.load(model_file)
        else:
            self.__model = None

        self.__mat_file = mat_file
        self.__outfile = outfile
        self.__best_geno = best_genotype
        self.__nodes = nodes
        self.__edges = edges
        self.__select_count = select_count
        self.__is_lower_better = is_lower_better
        self.__iter = iter
        self.__alpha = alpha
        self.__seed = seed
        self.best_alleles = []

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

    def model_predict_score(self, data):
        return self.__model.predict(np.array([data], dtype=np.int8))[0]

    def run(self):
        if self.__model is None:
            ssa = SelectSA(
                self.__select_count,
                self.__type_info,
                self.__best_geno,
                self.__mat_file,
                self.calc_score,
                self.__is_lower_better,
                iterate=self.__iter,
                alpha=self.__alpha,
                seed=self.__seed,
            )
        else:
            ssa = SelectSA(
                self.__select_count,
                self.__type_info,
                self.__best_geno,
                self.__mat_file,
                self.model_predict_score,
                self.__is_lower_better,
                iterate=self.__iter,
                alpha=self.__alpha,
                seed=self.__seed,
            )
        if ssa.run():
            for _ in range(len(ssa.data)):
                if ssa.data[_] == 1:
                    self.best_alleles.append(self.__allele_name[_])

            with open(self.__outfile, "w") as fout:
                fout.write("%s\n" % ("\n".join(self.best_alleles)))
            return True
        else:
            return False


class SelectNetSimple:
    def __init__(
        self,
        select_count,
        allele_name,
        outfile,
        best_genotype,
        nodes,
        edges,
        is_lower_better,
    ):
        self.__allele_name = allele_name
        self.__outfile = outfile
        self.__best_geno = best_genotype
        self.__nodes = nodes
        self.__edges = edges
        self.__select_count = select_count
        self.__is_lower_better = is_lower_better
        self.best_alleles = []

    def run(self):
        best_alleles = set()
        for idx in range(len(self.__best_geno)):
            if self.__best_geno[idx] == 1:
                best_alleles.add(self.__allele_name[idx])

        node_weights_db = {}
        for idx in range(len(self.__allele_name)):
            node_weights_db[self.__allele_name[idx]] = self.__nodes[idx]

        best_sel = []
        best_used_genes = set()
        for node in sorted(
            node_weights_db,
            key=lambda x: node_weights_db[x],
            reverse=False if self.__is_lower_better else True,
        ):
            if node not in best_alleles:
                continue
            gn, _ = node.split("-")
            if gn in best_used_genes:
                continue
            best_sel.append(node)
            best_used_genes.add(gn)
            if len(best_sel) >= self.__select_count:
                break

        if len(best_sel) >= self.__select_count:
            with open(self.__outfile, "w") as fout:
                fout.write("%s\n" % ("\n".join(best_sel)))
            return True
        return False
