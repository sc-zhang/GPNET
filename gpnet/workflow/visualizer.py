from gpnet.io.data_io import GraphLoader
from gpnet.io.message import Message
from pyecharts import options as opts
from pyecharts.charts import Graph


def plot_chart(nodes, edges, categories, out_html):
    g = (
        Graph(init_opts={"width": "1080px", "height": "1080px"})
        .add(
            "",
            nodes,
            edges,
            categories,
            repulsion=500,
            linestyle_opts=opts.LineStyleOpts(curve=0.2),
            label_opts=opts.LabelOpts(is_show=False)
        )
        .set_global_opts(
            legend_opts=opts.LegendOpts(is_show=False),
            title_opts=opts.TitleOpts(title="CWNet")
        )
        .render(out_html)
    )


def main(args):
    out_html = args.output
    in_data = args.input
    is_lower_better = args.lower

    Message.info("Loading data")
    dl = GraphLoader()
    dl.load_data(in_data, is_lower_better)

    Message.info("Plotting")
    plot_chart(dl.nodes, dl.edges, dl.categories, out_html)

    Message.info("Finished")
