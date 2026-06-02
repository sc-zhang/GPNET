from gpnet.algorithm.SA import SelectSA
import numpy as np
from gpnet.lib import fast_ridge


class SelectNet:
    def __init__(
        self,
        select_count,
        type_info,
        allele_name,
        modelfile,
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
        if modelfile is not None:
            self.__model = fast_ridge.GraphModel.load(modelfile)
        else:
            self.__model = None
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
                self.model_predict_score,
                self.__is_lower_better,
                iterate=self.__iter,
                alpha=self.__alpha,
                seed=self.__seed,
            )
        if ssa.run():
            for _ in range(len(ssa.geno)):
                if ssa.geno[_] == 1:
                    self.best_alleles.append(self.__allele_name[_])

            with open(self.__outfile, "w") as fout:
                fout.write("%s\n" % ("\n".join(self.best_alleles)))
            return True
        else:
            return False
