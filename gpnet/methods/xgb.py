import xgboost as xgb
from sklearn.model_selection import train_test_split
from gpnet.simulator.calc import calc_pheno
from gpnet.algorithm.SA import SA
from gpnet.io.data_io import DataLoader, DataSaver
from os import getpid


class XGB:
    def __init__(self, in_file, weight_file, out_file):
        self.__in_file = in_file
        self.__weight_file = weight_file
        self.__out_file = out_file

        self.__model = None

    def predict(self, data):
        return self.__model.predict(xgb.DMatrix([data]))[0]

    def __model_train(self, genotypes, phenotypes):
        params = {'learning_rate': 0.001,
                  'max_depth': 2,
                  'objective': 'reg:squarederror',
                  'gamma': 0,
                  'subsample': 0.7,
                  'colsample_bytree': 0.7,
                  'reg_alpha': 0.005,
                  'nthread': 6,
                  'eval_metric': ['logloss', 'rmse', 'mae'],
                  'eta': 0.3
                  }

        x_train, x_test, y_train, y_test = train_test_split(genotypes, phenotypes, test_size=.25)
        dtrain = xgb.DMatrix(x_train, label=y_train)
        dtest = xgb.DMatrix(x_test, label=y_test)

        res = xgb.cv(params, dtrain, num_boost_round=5000, metrics='rmse', early_stopping_rounds=25)
        best_nround = res.shape[0] - 1

        watchlist = [(dtrain, 'train'), (dtest, 'eval')]
        evals_result = {}

        self.__model = xgb.train(params,
                                 dtrain,
                                 num_boost_round=best_nround,
                                 evals=watchlist,
                                 evals_result=evals_result)

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

        print("\tPID:%d Starting XGBoost" % getpid())
        self.__model_train(genotypes, phenotypes)
        print("\tPID:%d Running SA" % getpid())
        sa = SA(type_info, self.predict, iterate=100, alpha=0.99)
        sa.run()
        best_sa_data = sa.data
        best_sa_pheno = calc_pheno(best_sa_data, single_weight, multi_weight)

        print("\tPID:%d Saving predict data" % getpid())
        ds = DataSaver(self.__out_file)
        ds.save_data(type_info, [best_sa_data], [best_sa_pheno])
