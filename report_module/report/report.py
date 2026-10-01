import os
from report_module.report.io import PopMatIO, PredictMatIO, Message, ReportIO


def get_predict_source(pop_mat_io: PopMatIO, predict_mat_io: PredictMatIO):
    source_db = {}
    sample_contain_allele_db = {}
    for allele in predict_mat_io.selected_alleles:
        if allele in pop_mat_io.allele_in_sample_idx_db:
            source_db[allele] = set()
            for sample_idx in sorted(pop_mat_io.allele_in_sample_idx_db[allele]):
                sample = pop_mat_io.idx_to_sample_db[sample_idx]
                source_db[allele].add(sample)
                if sample not in sample_contain_allele_db:
                    sample_contain_allele_db[sample] = set()
                sample_contain_allele_db[sample].add(allele)

    sample_score = {}
    for sample in sample_contain_allele_db:
        score = 0
        for allele in sample_contain_allele_db[sample]:
            if allele in predict_mat_io.node_weights:
                score += predict_mat_io.node_weights[allele]
        sample_score[sample] = [
            sample_contain_allele_db[sample],
            len(sample_contain_allele_db[sample]),
            score,
        ]

    return source_db, sample_score


def get_select_source(pop_mat_io: PopMatIO, select_list):
    source_db = {}
    sample_contain_allele_db = {}
    for allele in select_list:
        if allele in pop_mat_io.allele_in_sample_idx_db:
            source_db[allele] = set()
            for sample_idx in sorted(pop_mat_io.allele_in_sample_idx_db[allele]):
                sample = pop_mat_io.idx_to_sample_db[sample_idx]
                source_db[allele].add(sample)
                if sample not in sample_contain_allele_db:
                    sample_contain_allele_db[sample] = set()
                sample_contain_allele_db[sample].add(allele)
    sample_score = {
        sample: [
            sample_contain_allele_db[sample],
            len(sample_contain_allele_db[sample]),
        ]
        for sample in sample_contain_allele_db
    }

    return source_db, sample_score


def process(opts):
    input_file = opts.input
    pop_mat_file = opts.mat
    output_dir = opts.output

    os.makedirs(output_dir, exist_ok=True)

    Message.info("Checking input data")
    ft = ""
    with open(input_file, "r") as fin:
        for line in fin:
            if line[0] == "#":
                ft = "mat"
            else:
                ft = "sel"
            break

    Message.info("Input filetype is %s file" % ("predict" if ft == "mat" else "select"))

    Message.info("Loading population mat")
    pop_io = PopMatIO(pop_mat_file)
    pop_io.load()

    Message.info("Loading input data")
    if ft == "mat":
        pred_io = PredictMatIO(input_file)
        pred_io.load()
    else:
        select_alleles = []
        with open(input_file, "r") as fin:
            for line in fin:
                select_alleles.append(line.strip())

    Message.info("Generating report")
    if ft == "mat":
        source_db, sample_score = get_predict_source(pop_io, pred_io)
    else:
        source_db, sample_score = get_select_source(pop_io, select_alleles)

    Message.info("Writing report to file")
    allele_source_file = os.path.join(output_dir, "allele_source.tsv")
    rep_io = ReportIO(allele_source_file)
    rep_io.save_allele_source(source_db, pred_io.node_weights if ft == "mat" else None)

    sample_score_file = os.path.join(output_dir, "sample_score.tsv")
    rep_io = ReportIO(sample_score_file)
    rep_io.save_sample_score(sample_score, ft)

    Message.info("Finished")
