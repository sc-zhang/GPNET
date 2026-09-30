from gpnet.io.data_io import NetLoader
from gpnet.io.message import Message
from gpnet.selector.selector import SelectNet, SelectNetSimple


def main(args):
    in_data = args.input
    model_file = args.model
    mat_file = args.mat
    out_file = args.output
    select_cnt = args.count
    is_lower_better = args.lower
    method = args.method
    iter_cnt = args.iter
    iter_alpha = args.alpha
    seed = args.seed

    Message.info("Loading data: %s" % in_data)
    dl = NetLoader()
    dl.load_net(in_data)

    Message.info("Selecting with %d" % select_cnt)
    if select_cnt >= len(dl.type_info):
        Message.error(
            "Gene counts %d not greater than select count %d, Aborting..."
            % (len(dl.type_info), select_cnt)
        )
    else:
        if method == "SA":
            selector = SelectNet(
                select_cnt,
                dl.type_info,
                dl.allele_name,
                model_file,
                mat_file,
                out_file,
                dl.best_genotype,
                dl.nodes,
                dl.edges,
                is_lower_better,
                iter_cnt,
                iter_alpha,
                seed,
            )
            if selector.run():
                Message.info("Selected")
            else:
                Message.error("Unable select genes")
        elif method == "Simple":
            selector = SelectNetSimple(
                select_cnt,
                dl.allele_name,
                out_file,
                dl.best_genotype,
                dl.nodes,
                dl.edges,
                is_lower_better,
            )
            if selector.run():
                Message.info("Selected")
            else:
                Message.error("Unable select genes")
    Message.info("Finished")
