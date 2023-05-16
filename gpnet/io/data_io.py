from numpy import array


class DataSaver:
    def __init__(self, out_file):
        self.__out_file = out_file

    def save_data(self, type_info, genotypes, phenotypes, additional_info=None):
        # type_info is like:
        # [site1_type_count, site2_type_count, ..., siteN_type_count]
        with open(self.__out_file, 'w') as fout:
            fout.write("#TypeInfo\t%s\n" % ('\t'.join(map(str, type_info))))
            fout.write("#Sample\t")
            for _ in range(len(type_info)):
                for __ in range(type_info[_]):
                    fout.write("Site%d-Type%d\t" % (_ + 1, __ + 1))
            fout.write("Phenotype\n")
            for _ in range(len(genotypes)):
                fout.write("%d\t%s\t%s\n" % (_ + 1, '\t'.join(list(map(str, list(genotypes[_])))), str(phenotypes[_])))
            if additional_info:
                fout.write("%s\n" % ('\n'.join(additional_info)))

    def save_weight(self, single_weight, multi_weight):
        with open(self.__out_file, 'w') as fout:
            fout.write("## Single weight\n")
            fout.write("%s\n" % ('\t'.join(map(str, list(single_weight)))))
            fout.write("## Multi weight\n")
            for _ in sorted(multi_weight, reverse=True):
                for __ in multi_weight[_]:
                    fout.write("# %d\t%d\n" % (_, multi_weight[_][__]))
                    fout.write("%s\n" % ('\t'.join(map(str, list(__)))))


class DataLoader:
    def __init__(self):
        self.type_info = None
        self.genotypes = None
        self.phenotypes = None
        self.single_weight = None
        self.multi_weight = None

    def load_genotype(self, genotype_file):
        self.genotypes = []
        self.phenotypes = []
        with open(genotype_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                if line[0] == '#':
                    if line.startswith("#TypeInfo"):
                        self.type_info = list(map(int, data[1:]))
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
        with open(weight_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                if data[1] == 'Single':
                    is_single = True
                    continue
                elif data[1] == 'Multi':
                    is_single = False
                    continue
                if is_single:
                    self.single_weight = array(list(map(int, data)))
                else:
                    if line[0] == '#':
                        cur_pheno = int(data[2])
                        sn = int(data[1])
                        if sn not in self.multi_weight:
                            self.multi_weight[sn] = {}
                    else:
                        self.multi_weight[sn][tuple(sorted(list(map(int, data))))] = cur_pheno
