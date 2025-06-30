from gpnet.io.data_io import NetLoader
from gpnet.io.message import Message
from gpnet.selector.selector import SelectNet


def main(args):
    in_data = args.input
    out_file = args.output
    select_cnt = args.count
    is_lower_better = args.lower

    Message.info("Loading data")
    dl = NetLoader()
    dl.load_net(in_data)

    Message.info("Selecting")
    if select_cnt >= len(dl.type_info):
        Message.error("Gene counts %d not greater than select count %d, Aborting..." % (len(dl.type_info), select_cnt))
    else:
        selector = SelectNet(select_cnt, dl.type_info, dl.allele_name, out_file,
                             dl.best_genotype, dl.nodes, dl.edges, is_lower_better)
        if selector.run():
            Message.info("Selected")
        else:
            Message.error("Unable select genes")

    Message.info("Finished")
