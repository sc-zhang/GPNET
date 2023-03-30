from sklearn.linear_model import LinearRegression
from gpnet.simulator.calc import calc_pheno
from gpnet.algorithm.SA import SA
from gpnet.io.data_io import DataLoader, DataSaver
from os import getpid


class MLR:
    def __init__(self, in_file, weight_file, out_file):
        self.__in_file = in_file
        self.__weight_file = weight_file
        self.__out_file = out_file

        self.__model = None

    def predict(self, data):
        return self.__model.predict([data])[0]

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

        print("\tPID:%d Starting multi linear regression" % getpid())
        X_train = genotypes
        y_train = phenotypes
        lrg = LinearRegression()
        self.__model = lrg.fit(X_train, y_train)
        print("\tPID:%d Running SA" % getpid())
        sa = SA(type_info, self.predict, iterate=100, alpha=0.99)
        sa.run()
        best_sa_data = sa.data
        best_sa_pheno = calc_pheno(best_sa_data, single_weight, multi_weight)

        print("\tPID:%d Saving predict data" % getpid())
        ds = DataSaver(self.__out_file)
        ds.save_data(type_info, [best_sa_data], [best_sa_pheno])
