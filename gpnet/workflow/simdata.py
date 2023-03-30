from gpnet.io.data_io import DataSaver
from gpnet.simulator.simulator import Simulator
from os import getpid, getcwd, chdir, path, makedirs
from pathos.multiprocessing import Pool


def each_round(round_idx, site_cnt, max_type_cnt, sample_cnt):
    print("Round %d" % (round_idx + 1))
    print("\tPID:%d Generating genotypes" % getpid())
    simulator = Simulator(site_cnt, max_type_cnt, sample_cnt)
    simulator.sim_genotypes()
    genotypes = simulator.genotypes
    type_info = simulator.type_info

    print("\tPID:%d Generating genomic effects" % getpid())
    simulator.sim_site_effects()
    single_weight = simulator.single_weight
    multi_weight = simulator.multi_weight

    print("\tPID:%d Generating phenotypes" % getpid())
    simulator.sim_phenotypes()
    phenotypes = simulator.phenotypes

    print("\tPID:%d Saving data" % getpid())
    out_file = "Round%d.txt" % (round_idx + 1)
    ds = DataSaver(out_file)
    ds.save_data(type_info, genotypes, phenotypes)
    print("\tPID:%d Saving weight" % getpid())

    out_file = "Weight%d.txt" % (round_idx + 1)
    ds = DataSaver(out_file)
    ds.save_weight(single_weight, multi_weight)

    print("\tPID:%d Finished" % getpid())


def main(args):
    out_dir = args.outdir
    round_cnt = args.round
    site_cnt = args.sites
    max_type_cnt = args.max_types
    sample_cnt = args.samples
    thread = args.thread

    cur_dir = getcwd()
    if not path.exists(out_dir):
        makedirs(out_dir)
    chdir(out_dir)
    print("Simulating")
    pool = Pool(processes=thread)
    for _ in range(round_cnt):
        pool.apply_async(each_round, (_, site_cnt, max_type_cnt, sample_cnt, ))
    pool.close()
    pool.join()
    chdir(cur_dir)
    print("Finished")
