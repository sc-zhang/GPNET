from gpnet.methods.mlr import MLR
from gpnet.methods.cwnet import CWNET
from pathos.multiprocessing import Pool
from os import path, makedirs


def sub_process(idx, method, in_file, weight_file, out_file):
    print("Round %d" % (idx + 1))
    try:
        func = eval(method.upper())(in_file, weight_file, out_file)
        func.run()
    except NameError:
        print("No method named: %s" % method)


def load_file_list(in_file):
    file_list = []
    with open(in_file, 'r') as fin:
        for line in fin:
            file_list.append(line.strip())
    return file_list


def get_filename(full_file_path):
    return full_file_path.split('/')[-1].split('.')[0]


def main(args):
    in_file = args.input
    weight_file = args.weight
    output = args.output
    is_single = args.single
    method = args.method
    thread = args.thread

    print("Predicting")
    if is_single:
        sub_process(method, in_file, weight_file, output)
    else:
        in_file_list = load_file_list(in_file)
        weight_file_list = load_file_list(weight_file)
        if len(in_file_list) != len(weight_file_list):
            print("Input data and weight not match, please check it")
            exit(-1)
        if not path.exists(output):
            makedirs(output)
        pool = Pool(processes=thread)
        for i in range(len(in_file_list)):
            in_file = in_file_list[i]
            in_fn = get_filename(in_file)
            weight_file = weight_file_list[i]
            weight_fn = get_filename(weight_file)
            out_file = path.join(output, "%s-%s.txt" % (in_fn, weight_fn))
            pool.apply_async(sub_process, (i, method, in_file, weight_file, out_file,))
        pool.close()
        pool.join()
    print("Finished")
