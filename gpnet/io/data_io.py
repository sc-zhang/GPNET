from numpy import array
import numpy as np


class DataSaver:
    def __init__(self, out_file):
        self.__out_file = out_file

    def save_data(
        self, type_info, allele_name, genotypes, phenotypes, additional_info=None
    ):
        # type_info is like:
        # [site1_type_count, site2_type_count, ..., siteN_type_count]
        with open(self.__out_file, "w") as fout:
            fout.write("#TypeInfo\t%s\n" % ("\t".join(map(str, type_info))))
            fout.write("#Sample\t")
            for _ in allele_name:
                fout.write("%s\t" % _)
            fout.write("Phenotype\n")
            for _ in range(len(genotypes)):
                fout.write(
                    "%d\t%s\t%s\n"
                    % (
                        _ + 1,
                        "\t".join(list(map(str, list(genotypes[_])))),
                        str(phenotypes[_]),
                    )
                )
            if additional_info:
                fout.write("%s\n" % ("\n".join(additional_info)))

    def save_iter(self, type_info, allele_name, genotypes, score):
        # type_info is like:
        # [site1_type_count, site2_type_count, ..., siteN_type_count]
        with open(self.__out_file, "w") as fout:
            fout.write("#TypeInfo\t%s\n" % ("\t".join(map(str, type_info))))
            fout.write("#Iteration\t")
            for _ in allele_name:
                fout.write("%s\t" % _)
            fout.write("Score\n")
            for _ in range(len(score)):
                fout.write(
                    "%d\t%s\t%s\n"
                    % (
                        _ + 1,
                        "\t".join(list(map(str, list(genotypes[_])))),
                        str(score[_]),
                    )
                )

    def save_sim_data(self, type_info, genotypes, phenotypes):
        # type_info is like:
        # [site1_type_count, site2_type_count, ..., siteN_type_count]
        with open(self.__out_file, "w") as fout:
            fout.write("#TypeInfo\t%s\n" % ("\t".join(map(str, type_info))))
            fout.write("#Sample\t")
            for _ in range(len(type_info)):
                for __ in range(type_info[_]):
                    fout.write("Site%d-Type%d\t" % (_ + 1, __ + 1))
            fout.write("Phenotype\n")
            for _ in range(len(genotypes)):
                fout.write(
                    "%d\t%s\t%s\n"
                    % (
                        _ + 1,
                        "\t".join(list(map(str, list(genotypes[_])))),
                        str(phenotypes[_]),
                    )
                )

    def save_sim_weight(self, single_weight, multi_weight):
        with open(self.__out_file, "w") as fout:
            fout.write("## Single weight\n")
            fout.write("%s\n" % ("\t".join(map(str, list(single_weight)))))
            fout.write("## Multi weight\n")
            for _ in sorted(multi_weight, reverse=True):
                for __ in multi_weight[_]:
                    fout.write("# %d\t%d\n" % (_, multi_weight[_][__]))
                    fout.write("%s\n" % ("\t".join(map(str, list(__)))))


class DataLoader:
    def __init__(self):
        self.type_info = None
        self.allele_name = None
        self.genotypes = None
        self.phenotypes = None
        self.single_weight = None
        self.multi_weight = None

    def load_genotype(self, genotype_file):
        self.genotypes = []
        self.phenotypes = []
        with open(genotype_file, "r") as fin:
            for line in fin:
                data = line.strip().split()
                if line[0] == "#":
                    if line.startswith("#TypeInfo"):
                        self.type_info = list(map(int, data[1:]))
                    elif line.startswith("#Sample"):
                        self.allele_name = data[1:-1]
                    continue
                self.genotypes.append(list(map(int, data[1:-1])))
                self.phenotypes.append(float(data[-1]))
        self.genotypes = array(self.genotypes)
        self.phenotypes = array(self.phenotypes)

    def load_weight_data(self, weight_file):
        self.single_weight = []
        self.multi_weight = {}
        cur_pheno = 0
        is_single = False
        with open(weight_file, "r") as fin:
            for line in fin:
                data = line.strip().split()
                if data[1] == "Single":
                    is_single = True
                    continue
                elif data[1] == "Multi":
                    is_single = False
                    continue
                if is_single:
                    self.single_weight = array(list(map(int, data)))
                else:
                    if line[0] == "#":
                        cur_pheno = int(data[2])
                        sn = int(data[1])
                        if sn not in self.multi_weight:
                            self.multi_weight[sn] = {}
                    else:
                        self.multi_weight[sn][
                            tuple(sorted(list(map(int, data))))
                        ] = cur_pheno


class GraphLoader:
    def __init__(self):
        self.node_weights = {}
        self.edge_weights = []

    def load_data(
        self,
        in_file,
        is_lower_better,
        is_predicted_only,
    ):
        is_nodes = False
        is_edges = False
        predict_nodes = set()

        node_min = float("nan")
        node_max = float("nan")
        edge_min = float("nan")
        edge_max = float("nan")

        with open(in_file, "r") as fin:
            for line in fin:
                if line[0] == "#":
                    if line.startswith("#Sample"):
                        sample_list = line.strip().split()
                        continue
                    if line.startswith("# Nodes"):
                        is_nodes = True
                        continue
                    if line.startswith("# Edges"):
                        is_nodes = False
                        is_edges = True
                        continue
                    if is_nodes:
                        data = line.strip().split()
                        if len(data) < 3:
                            continue
                        smp = data[1]
                        if is_predicted_only and smp not in predict_nodes:
                            continue
                        val = float(data[-1])
                        if np.isnan(node_min) or val < node_min:
                            node_min = val
                        if np.isnan(node_max) or val > node_max:
                            node_max = val
                        self.node_weights[smp] = val
                    elif is_edges:
                        data = line.strip().split()
                        if len(data) < 4:
                            continue
                        src = data[1]
                        tgt = data[2]
                        val = float(data[-1])
                        if is_predicted_only and (
                            src not in predict_nodes or tgt not in predict_nodes
                        ):
                            continue
                        if np.isnan(edge_min) or val < edge_min:
                            edge_min = val
                        if np.isnan(edge_max) or val > edge_max:
                            edge_max = val
                        self.edge_weights.append((src, tgt, val))
                else:
                    data = line.strip().split()
                    for _ in range(1, len(data) - 1):
                        if data[_] == "1":
                            predict_nodes.add(sample_list[_])


class GraphSaver:
    def __init__(self, out_pre):
        self.__out_pre = out_pre

    def save_graph(self, graph):
        with open(self.__out_pre + ".nodes.csv", "w") as fout:
            fout.write("Nodes,Weights\n")
            for _ in graph.node_weights:
                fout.write("%s,%f\n" % (_, graph.node_weights[_]))

        with open(self.__out_pre + ".edges.csv", "w") as fout:
            fout.write("Source,Target,Weight\n")
            for src, tgt, val in graph.edge_weights:
                fout.write("%s,%s,%f\n" % (src, tgt, val))


class NetLoader:
    def __init__(self):
        self.type_info = None
        self.allele_name = None
        self.best_genotype = None
        self.nodes = None
        self.edges = None

    def load_net(self, in_file):
        is_node = False
        is_edge = False
        with open(in_file, "r") as fin:
            for line in fin:
                if line[0] == "#":
                    data = line.strip().split()
                    if line.startswith("#TypeInfo"):
                        self.type_info = list(map(int, data[1:]))
                        continue
                    if line.startswith("#Sample"):
                        self.allele_name = data[1:-1]
                        allele_idx = {
                            self.allele_name[_]: _ for _ in range(len(self.allele_name))
                        }
                        self.nodes = [0 for _ in range(len(self.allele_name))]
                        self.edges = [
                            [0 for _ in range(len(self.allele_name))]
                            for __ in range(len(self.allele_name))
                        ]
                        continue
                    if line.startswith("# Nodes"):
                        is_node = True
                        continue
                    if line.startswith("# Edges"):
                        is_node = False
                        is_edge = True
                        continue
                    if is_node:
                        if len(data) < 3:
                            continue
                        gid = data[1]
                        val = float(data[2])
                        self.nodes[allele_idx[gid]] = val
                    if is_edge:
                        if len(data) < 4:
                            continue
                        g1, g2, val = data[1:]
                        val = float(val)
                        gidx1 = allele_idx[g1]
                        gidx2 = allele_idx[g2]
                        self.edges[gidx1][gidx2] = val
                        self.edges[gidx2][gidx1] = val
                else:
                    data = line.strip().split("\t")
                    self.best_genotype = array(list(map(int, data[1:-1])))
