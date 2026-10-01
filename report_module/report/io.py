import time


class Message:
    def __init__(self):
        pass

    @staticmethod
    def info(info):
        print(
            "\033[32m%s\033[0m %s"
            % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
        )

    @staticmethod
    def warn(info):
        print(
            "\033[33m%s\033[0m %s"
            % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
        )

    @staticmethod
    def error(info):
        print(
            "\033[31m%s\033[0m %s"
            % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
        )


class PopMatIO:
    def __init__(self, pop_mat_file):
        self.__pop_mat_file = pop_mat_file
        self.idx_to_allele_db = {}
        self.allele_in_sample_idx_db = {}
        self.idx_to_sample_db = {}

    def load(self):
        sample_idx = 0
        with open(self.__pop_mat_file, "r") as fin:
            for line in fin:
                if line[0] == "#":
                    if line.startswith("#Sample"):
                        data = line.strip().split("\t")
                        for idx in range(1, len(data) - 1):
                            allele_name = data[idx]
                            self.idx_to_allele_db[idx] = allele_name
                            self.allele_in_sample_idx_db[allele_name] = set()
                else:
                    data = line.strip().split("\t")
                    sample_name = data[0]
                    self.idx_to_sample_db[sample_idx] = sample_name

                    for idx in range(1, len(data) - 1):
                        if data[idx] == "1":
                            allele_name = self.idx_to_allele_db[idx]
                            self.allele_in_sample_idx_db[allele_name].add(sample_idx)

                    sample_idx += 1


class PredictMatIO:
    def __init__(self, predict_mat_file):
        self.__predict_mat_file = predict_mat_file
        self.node_weights = {}
        self.selected_alleles = []

    def load(self):
        is_nodes = False
        with open(self.__predict_mat_file, "r") as fin:
            for line in fin:
                if line[0] == "#":
                    if line.startswith("#Sample"):
                        data = line.strip().split("\t")
                        allele_list = data[1:-1]
                        continue
                    if line.startswith("# Nodes"):
                        is_nodes = True
                        continue
                    if line.startswith("# Edges"):
                        is_nodes = False
                        continue
                    if is_nodes:
                        data = line.strip().split()
                        if len(data) < 3:
                            continue
                        allele_name = data[1]
                        weight = float(data[-1])
                        self.node_weights[allele_name] = weight
                else:
                    data = line.strip().split()
                    for idx in range(1, len(data) - 1):
                        if data[idx] == "1":
                            self.selected_alleles.append(allele_list[idx - 1])


class ReportIO:
    def __init__(self, output_file):
        self.__output_file = output_file

    def save_allele_source(self, data_db, weight_db=None):
        if weight_db:
            with open(self.__output_file, "w") as fout:
                fout.write("#AlleleName\tWeight\tSourceSamples\n")
                for allele_name in sorted(
                    data_db, key=lambda x: weight_db[x], reverse=True
                ):
                    fout.write(
                        "%s\t%f\t%s\n"
                        % (
                            allele_name,
                            weight_db[allele_name],
                            ",".join(sorted(data_db[allele_name])),
                        )
                    )
        else:
            with open(self.__output_file, "w") as fout:
                fout.write("#AlleleName\tSourceSamples\n")
                for allele_name in sorted(data_db):
                    fout.write(
                        "%s\t%s\n"
                        % (allele_name, ",".join(sorted(data_db[allele_name])))
                    )

    def save_sample_score(self, data_db, input_type="mat"):
        with open(self.__output_file, "w") as fout:
            if input_type == "mat":
                fout.write("#Sample\tContainAlleleCount\tScore\tAlleles\n")
            else:
                fout.write("#Sample\tContainAlleleCount\tAlleles\n")

            for sample in sorted(data_db, key=lambda x: data_db[x][-1], reverse=True):
                if input_type == "mat":
                    fout.write(
                        "%s\t%d\t%f\t%s\n"
                        % (
                            sample,
                            data_db[sample][1],
                            data_db[sample][2],
                            ",".join(sorted(data_db[sample][0])),
                        )
                    )
                else:
                    fout.write(
                        "%s\t%d\t%s\n"
                        % (
                            sample,
                            data_db[sample][1],
                            ",".join(sorted(data_db[sample][0])),
                        )
                    )
